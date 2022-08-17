#!/usr/bin/env python
# Copyright 2003, by Kristian Rother. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# PDBList.py
#
# A tool for tracking changes in the PDB Protein Structure Database.
#
# (c) 2003 Kristian Rother
# This work was supported by the German Ministry of Education
# and Research (BMBF). Project http://www.bcbio.de
#
# Contact the author
#    homepage : http://www.rubor.de/bioinf
#    email    : krother@genesilico.pl
#
#
# (c) 2016 Wiktoria Karwicka & Jacek Smietanski
#   - updated and Python 3.x compatible code
#   - new options to enable download PDBx/mmCif, PDBML and mmtf formatted
#       files as well as large PDB bundles
#   - unit tests for the module
#
# Contact the corresponding author
#   homepage : http://jaceksmietanski.net
#   email    : jacek.smietanski@ii.uj.edu.pl
#
# It may be distributed freely with respect to the original authors.
# Any maintainer of the Biopython code may change this notice
# when appropriate.

"""Access the PDB over the internet (e.g. to download structures)."""

from __future__ import annotations

import contextlib
import dataclasses
import functools
import gzip
import os
import re
import shutil
import socket
import sys
import time
from urllib.error import URLError
from urllib.parse import urljoin
from urllib.parse import urlsplit
from urllib.request import urlcleanup
from urllib.request import urlopen
from urllib.request import urlretrieve

from Bio.PDB.PDBExceptions import PDBException


@dataclasses.dataclass
class PDBServer:
    """Dataclass to store parameters for connecting to PDB server.

    param pdb_dir_url: URL to the pdb directory.
    type pdb_dir_url: string
    """

    pdb_dir_url: str

    @functools.cached_property
    def entries(self) -> list[str]:
        """List of PDB codes based on entries.idx."""
        url = urljoin(self.pdb_dir_url, "derived_data/index/entries.idx")
        return re.findall(r"(\w{4})\t.*", read_url(url).decode())

    @functools.cached_property
    def latests(self) -> tuple[list[str], list[str], list[str]]:
        """Return three lists of the newest weekly files (added,mod,obsolete)."""
        latest_dir_url = urljoin(self.pdb_dir_url, "data/status/latest/")
        latests = {}
        latest_regex = re.compile(r"(\w{4})\n.*")
        for latest in ["added", "modified", "obsolete"]:
            latest_url = urljoin(latest_dir_url, f"{latest}.pdb")
            latests[latest] = latest_regex.findall(read_url(latest_url).decode())
        return (latests["added"], latests["modified"], latests["obsolete"])

    @functools.cached_property
    def obsoletes(self) -> list[str]:
        """Return a list of all obsolete entries ever in the PDB.

        Returns a list of all obsolete pdb codes that have ever been
        in the PDB.

        Gets and parses the file from the PDB server in the format
        (the first pdb_code column is the one used). The file looks
        like this::

             LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS
            OBSLTE    31-JUL-94 116L     216L
            ...
            OBSLTE    29-JAN-96 1HFT     2HFT
            OBSLTE    21-SEP-06 1HFV     2J5X
            OBSLTE    21-NOV-03 1HG6
            OBSLTE    18-JUL-84 1HHB     2HHB 3HHB
            OBSLTE    08-NOV-96 1HID     2HID
            OBSLTE    01-APR-97 1HIU     2HIU
            OBSLTE    14-JAN-04 1HKE     1UUZ
            ...

        """
        url = urljoin(self.pdb_dir_url, "data/status/obsolete.dat")
        return re.findall(
            r"OBSLTE\s+\d{2}-\w{3}-\d{2}\s+(\w{4})", read_url(url).decode()
        )

    @functools.cached_property
    def sequences(self):
        """Retrieve and save a (big) file containing all the sequences of PDB entries."""
        url = urljoin(self.pdb_dir_url, "derived_data/pdb_seqres.txt")
        return read_url(url)


@functools.lru_cache(maxsize=None)
def read_url(url: str):
    """Read resource from URL."""
    with contextlib.closing(urlopen(url)) as handle:
        return handle.read()


SERVERS = [
    PDBServer("https://files.rcsb.org/pub/pdb/"),
    PDBServer("https://s3.rcsb.org/pub/pdb/"),
    PDBServer("https://ftp.ebi.ac.uk/pub/databases/pdb/"),
    PDBServer("https://ftp.pdbj.org/pub/pdb/"),
]


class PDBList:
    """Quick access to the structure lists on the PDB or its mirrors.

    This class provides quick access to the structure lists on the
    PDB server or its mirrors. The structure lists contain
    four-letter PDB codes, indicating that structures are
    new, have been modified or are obsolete. The lists are released
    on a weekly basis.

    It also provides a function to retrieve PDB files from the server.
    To use it properly, prepare a directory /pdb or the like,
    where PDB files are stored.

    All available file formats (PDB, PDBx/mmCif, PDBML, mmtf) are supported.
    Please note that large structures (containing >62 chains
    and/or 99999 ATOM lines) are no longer stored as a single PDB file
    and by default (when PDB format selected) are not downloaded.

    Large structures can be downloaded in other formats, including PDBx/mmCif
    or as a .tar file (a collection of PDB-like formatted files for a given
    structure).

    If you want to use this module from inside a proxy, add
    the proxy variable to your environment, e.g. in Unix:
    export HTTP_PROXY='http://realproxy.charite.de:888'
    (This can also be added to ~/.bashrc)
    """

    PDB_REF = """
    The Protein Data Bank: a computer-based archival file for macromolecular structures.
    F.C.Bernstein, T.F.Koetzle, G.J.B.Williams, E.F.Meyer Jr, M.D.Brice, J.R.Rodgers, O.Kennard, T.Shimanouchi, M.Tasumi
    J. Mol. Biol. 112 pp. 535-542 (1977)
    http://www.pdb.org/.
    """

    def __init__(
        self, server: str | None = None, pdb=None, obsolete_pdb=None, verbose=True
    ):
        """Initialize the class with the default server or a custom one.

        Argument pdb is the local path to use, defaulting to the current
        directory at the moment of initialisation.
        """
        self.pdb_server = PDBServer(server) if server else get_fastest_server()

        if pdb:
            self.local_pdb = pdb  # local pdb file tree
        else:
            self.local_pdb = os.getcwd()

        # enable or disable verbose
        self._verbose = verbose

        # local file tree for obsolete pdb files
        if obsolete_pdb:
            self.obsolete_pdb = obsolete_pdb
        else:
            self.obsolete_pdb = os.path.join(self.local_pdb, "obsolete")
            if not os.access(self.obsolete_pdb, os.F_OK):
                os.makedirs(self.obsolete_pdb)

        # variable for command-line option
        self.flat_tree = False

    @staticmethod
    def _print_default_format_warning(file_format):
        """Print a warning to stdout (PRIVATE).

        Temporary warning (similar to a deprecation warning) that files
        are being downloaded in mmCIF.
        """
        if file_format is None:
            sys.stderr.write(
                "WARNING: The default download format has changed from PDB to PDBx/mmCif\n"
            )
            return "mmCif"
        return file_format

    def retrieve_pdb_file(
        self, pdb_code, obsolete=False, pdir=None, file_format=None, overwrite=False
    ):
        """Fetch PDB structure file from PDB server, and store it locally.

        The PDB structure's file name is returned as a single string.
        If obsolete ``==`` True, the file will be saved in a special file tree.

        NOTE. The default download format has changed from PDB to PDBx/mmCif

        :param pdb_code: 4-symbols structure Id from PDB (e.g. 3J92).
        :type pdb_code: string

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PDBML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure)

        :type file_format: string

        :param overwrite: if set to True, existing structure files will be overwritten. Default: False
        :type overwrite: bool

        :param obsolete:
            Has a meaning only for obsolete structures. If True, download the obsolete structure
            to 'obsolete' folder, otherwise download won't be performed.
            This option doesn't work for mmtf format as obsoleted structures aren't stored in mmtf.
            Also doesn't have meaning when parameter pdir is specified.
            Note: make sure that you are about to download the really obsolete structure.
            Trying to download non-obsolete structure into obsolete folder will not work
            and you face the "structure doesn't exists" error.
            Default: False

        :type obsolete: bool

        :param pdir: put the file in this directory (default: create a PDB-style directory tree)
        :type pdir: string

        :return: filename
        :rtype: string
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)

        # Get the compressed PDB structure
        code = pdb_code.lower()
        archive = {
            "pdb": "pdb%s.ent.gz",
            "mmCif": "%s.cif.gz",
            "xml": "%s.xml.gz",
            "mmtf": "%s",
            "bundle": "%s-pdb-bundle.tar.gz",
        }
        archive_fn = archive[file_format] % code

        if file_format not in archive.keys():
            raise (
                "Specified file_format %s doesn't exists or is not supported. Maybe a "
                "typo. Please, use one of the following: mmCif, pdb, xml, mmtf, bundle"
                % file_format
            )

        if file_format in ("pdb", "mmCif", "xml"):
            pdb_dir = "divided" if not obsolete else "obsolete"
            file_type = (
                "pdb"
                if file_format == "pdb"
                else "mmCIF"
                if file_format == "mmCif"
                else "XML"
            )
            url = urljoin(
                self.pdb_server.pdb_dir_url,
                f"data/structures/{pdb_dir}/{file_type}/{code[1:3]}/{archive_fn}",
            )
        elif file_format == "bundle":
            url = urljoin(
                self.pdb_server.pdb_dir_url,
                f"compatible/pdb_bundle/{code[1:3]}/{code}/{archive_fn}",
            )
        else:
            url = f"http://mmtf.rcsb.org/v1.0/full/{code}"

        # Where does the final PDB file get saved?
        if pdir is None:
            path = self.local_pdb if not obsolete else self.obsolete_pdb
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = os.path.join(path, code[1:3])
        else:  # Put in specified directory
            path = pdir
        if not os.access(path, os.F_OK):
            os.makedirs(path)
        filename = os.path.join(path, archive_fn)
        final = {
            "pdb": "pdb%s.ent",
            "mmCif": "%s.cif",
            "xml": "%s.xml",
            "mmtf": "%s.mmtf",
            "bundle": "%s-pdb-bundle.tar",
        }
        final_file = os.path.join(path, final[file_format] % code)

        # Skip download if the file already exists
        if not overwrite:
            if os.path.exists(final_file):
                if self._verbose:
                    print(f"Structure exists: '{final_file}' ")
                return final_file

        # Retrieve the file(s)
        if self._verbose:
            print(f"Downloading PDB structure '{pdb_code}'...")
        try:
            urlcleanup()
            urlretrieve(url, filename)
        except OSError:
            print("Desired structure doesn't exist")
        else:
            with gzip.open(filename, "rb") as gz:
                with open(final_file, "wb") as out:
                    out.writelines(gz)
            os.remove(filename)
        return final_file

    def update_pdb(self, file_format=None, with_assemblies=False):
        """Update your local copy of the PDB files.

        I guess this is the 'most wanted' function from this module.
        It gets the weekly lists of new and modified pdb entries and
        automatically downloads the according PDB files.
        You can call this module as a weekly cron job.
        """
        assert os.path.isdir(self.local_pdb)
        assert os.path.isdir(self.obsolete_pdb)

        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)

        new, modified, obsolete = self.pdb_server.latests

        for pdb_code in new + modified:
            try:
                self.retrieve_pdb_file(pdb_code, file_format=file_format)
                if with_assemblies:
                    self.download_assemblies(
                        pdb_code,
                        file_format=file_format,
                        overwrite=True,
                    )
            except Exception as err:
                print(f"error {pdb_code}: {err}\n")
                # you can insert here some more log notes that
                # something has gone wrong.

        # Move the obsolete files to a special folder
        # NOTE: This should be updated to handle multiple file types and
        # assemblies. As of now, it only looks for PDB-formatted files.
        # Using pathlib will be probably the best approach here, to build
        # and index of which files we have efficiently (or glob them).
        for pdb_code in obsolete:
            if self.flat_tree:
                old_file = os.path.join(self.local_pdb, f"pdb{pdb_code}.ent")
                new_dir = self.obsolete_pdb
            else:
                old_file = os.path.join(
                    self.local_pdb, pdb_code[1:3], f"pdb{pdb_code}.ent"
                )
                new_dir = os.path.join(self.obsolete_pdb, pdb_code[1:3])
            new_file = os.path.join(new_dir, f"pdb{pdb_code}.ent")
            if os.path.isfile(old_file):
                if not os.path.isdir(new_dir):
                    os.mkdir(new_dir)
                try:
                    shutil.move(old_file, new_file)
                except Exception:
                    print(f"Could not move {old_file} to obsolete folder")
            elif os.path.isfile(new_file):
                if self._verbose:
                    print(f"Obsolete file {old_file} already moved")
            else:
                if self._verbose:
                    print(f"Obsolete file {old_file} is missing")

    def download_pdb_files(
        self, pdb_codes, obsolete=False, pdir=None, file_format=None, overwrite=False
    ):
        """Fetch set of PDB structure files from the PDB server and stores them locally.

        The PDB structure's file name is returned as a single string.
        If obsolete ``==`` True, the files will be saved in a special file tree.

        :param pdb_codes: a list of 4-symbols structure Ids from PDB
        :type pdb_codes: list of strings

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PMDML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure)

        :param overwrite: if set to True, existing structure files will be overwritten. Default: False
        :type overwrite: bool

        :param obsolete:
            Has a meaning only for obsolete structures.
            If True, download the obsolete structure
            to 'obsolete' folder, otherwise download won't be performed.
            This option doesn't work for mmtf format as obsoleted structures are not available as mmtf.
            (default: False)

        :type obsolete: bool

        :param pdir: put the file in this directory (default: create a PDB-style directory tree)
        :type pdir: string

        :return: filenames
        :rtype: string
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        for pdb_code in pdb_codes:
            self.retrieve_pdb_file(
                pdb_code,
                obsolete=obsolete,
                pdir=pdir,
                file_format=file_format,
                overwrite=overwrite,
            )

    def retrieve_assembly_file(
        self, pdb_code, assembly_num, pdir=None, file_format=None, overwrite=False
    ):
        """Fetch one or more assembly structures associated with a PDB entry.

        Unless noted below, parameters are described in ``retrieve_pdb_file``.

        :type  assembly_num: int
        :param assembly_num: assembly number to download.

        :rtype : str
        :return: file name of the downloaded assembly file.
        """
        archive = {
            "pdb": "%s.pdb%s.gz",
            "mmcif": "%s-assembly%s.cif.gz",
        }

        file_format = self._print_default_format_warning(file_format)
        file_format = file_format.lower()  # we should standardize this.
        if file_format not in archive:
            raise (
                "Specified file_format '%s' is not supported. Use one of the "
                "following: 'mmcif' or 'pdb'." % file_format
            )

        # Get the compressed assembly structure name
        archive_fn = archive[file_format] % (pdb_code.lower(), int(assembly_num))

        if file_format == "mmcif":
            url = urljoin(
                self.pdb_server.pdb_dir_url, f"data/assemblies/mmCIF/all/{archive_fn}"
            )
        elif file_format == "pdb":
            url = urljoin(
                self.pdb_server.pdb_dir_url, f"data/biounit/PDB/all/{archive_fn}"
            )
        else:  # better safe than sorry
            raise ValueError("file_format '%s' not supported: %s" % file_format)

        # Where will the file be saved?
        if pdir is None:
            path = self.local_pdb
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = os.path.join(path, pdb_code[1:3])
        else:  # Put in specified directory
            path = pdir
        if not os.access(path, os.F_OK):
            os.makedirs(path)

        assembly_gz_file = os.path.join(path, archive_fn)
        assembly_final_file = os.path.join(path, archive_fn[:-3])  # no .gz

        # Skip download if the file already exists
        if not overwrite:
            if os.path.exists(assembly_final_file):
                if self._verbose:
                    print(f"Structure exists: '{assembly_final_file}' ")
                return assembly_final_file

        # Otherwise,retrieve the file(s)
        if self._verbose:
            print(
                f"Downloading assembly ({assembly_num}) for PDB entry "
                f"'{pdb_code}'..."
            )
        urlcleanup()
        urlretrieve(url, assembly_gz_file)
        with gzip.open(assembly_gz_file, "rb") as gz:
            with open(assembly_final_file, "wb") as out:
                out.writelines(gz)
        os.remove(assembly_gz_file)
        return assembly_final_file

    def download_assemblies(
        self, pdb_code: str, file_format: str, overwrite: bool, pdir: str | None = None
    ):
        """Download all available assemblies for a given pdb_code."""
        max_index = 20  # safeguard
        assembly_index = 0
        next_assembly = True
        while next_assembly and assembly_index < max_index:
            assembly_index += 1
            try:
                self.retrieve_assembly_file(
                    pdb_code=pdb_code,
                    assembly_num=assembly_index,
                    file_format=file_format,
                    overwrite=overwrite,
                    pdir=pdir,
                )
            except URLError:
                next_assembly = False

    def download_all_assemblies(self, file_format=None):
        """Retrieve all biological assemblies not in the local PDB copy.

        :type  listfile: str, optional
        :param listfile: file name to which all assembly codes will be written

        :type  file_format: str, optional
        :param file_format: format in which to download the entries. Available
            options are "mmCif" or "pdb". Defaults to mmCif.
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        for pdb_code in self.pdb_server.entries:
            self.download_assemblies(pdb_code, file_format)

    def download_entire_pdb(self, listfile=None, file_format=None):
        """Retrieve all PDB entries not present in the local PDB copy.

        :param listfile: filename to which all PDB codes will be written (optional)

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PMDML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure)

        NOTE. The default download format has changed from PDB to PDBx/mmCif
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        for pdb_code in self.pdb_server.entries:
            self.retrieve_pdb_file(pdb_code, file_format=file_format)
        # Write the list
        if listfile:
            with open(listfile, "w") as outfile:
                outfile.write("\n".join(self.pdb_server.entries))

    def download_obsolete_entries(self, listfile=None, file_format=None):
        """Retrieve all obsolete PDB entries not present in local obsolete PDB copy.

        :param listfile: filename to which all PDB codes will be written (optional)

        :param file_format: file format. Available options:
            "mmCif" (default, PDBx/mmCif file),
            "pdb" (format PDB),
            "xml" (PMDML/XML format),

        NOTE. The default download format has changed from PDB to PDBx/mmCif
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        for pdb_code in self.pdb_server.obsoletes:
            self.retrieve_pdb_file(pdb_code, obsolete=True, file_format=file_format)

        # Write the list
        if listfile:
            with open(listfile, "w") as outfile:
                outfile.writelines("\n".join(self.pdb_server.obsoletes))

    def get_seqres_file(self, savefile="pdb_seqres.txt"):
        """Retrieve and save a (big) file containing all the sequences of PDB entries."""
        if self._verbose:
            print("Retrieving sequence file (takes over 110 MB).")
        with open(savefile, "wb+") as steam:
            steam.write(self.pdb_server.sequences)


@functools.lru_cache(maxsize=None)
def get_fastest_server(servers: list[PDBServer] | None = None) -> PDBServer:
    """Fastest server among the available ones."""
    servers = servers or SERVERS
    fastest_server = None
    fastest_timing = None

    for server in servers:
        url_parts = urlsplit(server.pdb_dir_url)

        server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        server_socket.settimeout(1)

        time_before_connect = time.time()
        try:
            server_socket.connect((url_parts.hostname, url_parts.port or 443))
        except Exception:
            continue
        finally:
            server_socket.close()
        timing = time.time() - time_before_connect

        if timing and (not fastest_timing or timing < fastest_timing):
            fastest_server = server
            fastest_timing = timing

    if not fastest_server:
        raise PDBException("Unable to connect to any PDB servers.")

    return fastest_server


if __name__ == "__main__":

    doc = """PDBList.py
    (c) Kristian Rother 2003, Wiktoria Karwicka & Jacek Smietanski 2016
    Contributed to Biopython

    Usage::

        PDBList.py update <pdb_path> [options]   - write weekly PDB updates to
                                                   local pdb tree.
        PDBList.py all    <pdb_path> [options]   - write all PDB entries to
                                                   local pdb tree.
        PDBList.py obsol  <pdb_path> [options]   - write all obsolete PDB
                                                   entries to local pdb tree.
        PDBList.py assemb <pdb_path> [options]   - write all assemblies for each
                                                   PDB entry to local pdb tree.
        PDBList.py <PDB-ID> <pdb_path> [options] - retrieve single structure
        PDBList.py (<PDB-ID1>,<PDB-ID2>,...) <pdb_path> [options] - retrieve a set
                                                   of structures

    Options:
     -d       A single directory will be used as <pdb_path>, not a tree.
     -o       Overwrite existing structure files.
     -pdb     Downloads structures in PDB format
     -xml     Downloads structures in PDBML (XML) format
     -mmtf    Downloads structures in mmtf format
     -with-assemblies    Downloads assemblies along with regular entries.

    Maximum one format can be specified simultaneously (if more selected, only
    the last will be considered). By default (no format specified) structures are
    downloaded as PDBx/mmCif files.
    """
    print(doc)

    file_format = "mmCif"
    overwrite = False
    with_assemblies = False

    if len(sys.argv) > 2:
        pdb_path = sys.argv[2]
        pl = PDBList(pdb=pdb_path)
        if len(sys.argv) > 3:
            for option in sys.argv[3:]:
                if option == "-d":
                    pl.flat_tree = True
                elif option == "-o":
                    overwrite = True
                elif option in ("-pdb", "-xml", "-mmtf"):
                    file_format = option[1:]
                # Allow for download of assemblies alongside ASU
                elif option == "-with-assemblies":
                    with_assemblies = True

    else:
        pdb_path = os.getcwd()
        pl = PDBList()
        pl.flat_tree = True

    if len(sys.argv) > 1:
        if sys.argv[1] == "update":
            # update PDB
            print("updating local PDB at " + pdb_path)
            pl.update_pdb(file_format=file_format, with_assemblies=with_assemblies)

        elif sys.argv[1] == "all":
            # get the entire PDB
            pl.download_entire_pdb(file_format=file_format)
            if with_assemblies:
                # get all assembly structures
                pl.download_all_assemblies(file_format=file_format)

        elif sys.argv[1] == "obsol":
            # get all obsolete entries
            pl.download_obsolete_entries(pdb_path, file_format=file_format)

        elif sys.argv[1] == "assemb":
            # get all assembly structures
            pl.download_all_assemblies(file_format=file_format)

        elif len(sys.argv[1]) == 4 and sys.argv[1][0].isdigit():
            pdb_code = sys.argv[1]
            # get single PDB entry
            pl.retrieve_pdb_file(
                pdb_code, pdir=pdb_path, file_format=file_format, overwrite=overwrite
            )
            if with_assemblies:
                # PDB Code might have more than one assembly.
                pl.download_assemblies(
                    pdb_code=pdb_code,
                    file_format=file_format,
                    overwrite=overwrite,
                )

        elif sys.argv[1][0] == "(":
            # get a set of PDB entries
            pdb_ids = re.findall("[0-9A-Za-z]{4}", sys.argv[1])
            for pdb_id in pdb_ids:
                pl.retrieve_pdb_file(
                    pdb_id, pdir=pdb_path, file_format=file_format, overwrite=overwrite
                )
                if with_assemblies:
                    # PDB Code might have more than one assembly.
                    assemblies = pl.download_assemblies(
                        pdb_code=pdb_id,
                        file_format=file_format,
                        overwrite=overwrite,
                    )
