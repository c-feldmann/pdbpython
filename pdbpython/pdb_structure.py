import os
import warnings
from typing import *
import urllib3
import time
from collections import defaultdict

AMINO_ACID2LETTER = ["ala", "cys", "asp", "glu", "phe", "gly", "his", "ile", "lys", "leu", "met", "asn", "pyl", "pro",
                     "gln", "arg", "ser", "thr", "sec", "val", "trp", "tyr"]


class MissingResidueError(Exception):
    pass


def query_website(request: str, max_trials: int = 10):
    """A website specified in request is queried repeatedly until success or max_trails is exceeded.

    Parameters
    ----------
    request: str
        website plus query.
    max_trials
        maximum of trials.
    Returns
    -------
        encoded string.
    Raises
    ______
    ValueError
        When the return code is equal to 404 the URL itself might be faulty. No repeats are executed.
    ConnectionError
        If after `max_trials` no successful connection was established.
    """
    response = None
    for trial in range(max_trials):
        with urllib3.PoolManager() as pool:
            response = pool.request('GET', request)
            if response.status == 200:
                break
            elif response.status == 404:
                raise ValueError(f"Invalid url: {request}")
        time.sleep(trial * 10)
    else:
        raise ConnectionError(response.status)
    return response.data


class PDBAtom:
    """ Contains information of an 'ATOM' or 'HETATM' line in a PDB file.

    Direct creation of such an object is possible but not encouraged. Best practise is to create these objects from a
    line in a PDB file with:
        atom = PDBAtom.from_pdb_line(line)

    Attributes
    __________
    line: str
        reconstructed line from atom properties.
    x: float
        x-coordinate
    y: float
        y-coordinate
    z: float
        z-coordinate
    """
    def __init__(self,
                 is_het: bool,
                 atom_id: int,
                 atom_label: str,
                 alt_pos: str,
                 res_name: str,
                 chain: str,
                 res_id: int,
                 icode: str,
                 coordinates: List[float],
                 occupancy: Optional[float],
                 temperature_factor: Optional[float],
                 segment_id: str,
                 element: str,
                 charge: Optional[str],
                 ):
        """ Initializes the atom.

        For the meaning of parameters please refer to:
            www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

        Parameters
        ----------
         is_het: bool
         atom_id: int
         atom_label: str
         alt_pos: str
         res_name: str
         chain: str
         res_id: int
         coordinates: List[float]
         occupancy: Optional[float]
         temperature_factor: Optional[float]
         segment_id: str
         element: str
         charge: Optional[str]
        """
        self.is_het = is_het
        self.atom_id = atom_id
        self.atom_label = atom_label
        self.alt_pos = alt_pos
        self.res_name = res_name
        self.chain = chain
        self.res_id = res_id
        self.icode = icode
        self.coordinates = coordinates
        self._occupancy = occupancy
        self._temperature_factor = temperature_factor
        self.segment_id = segment_id
        self.element = element
        self.charge = charge

        self._line = None

    @property
    def line(self):
        return self._reconstructed_line

    @property
    def x(self) -> float:
        return self.coordinates[0]

    @property
    def y(self) -> float:
        return self.coordinates[1]

    @property
    def z(self) -> float:
        return self.coordinates[2]

    @property
    def occupancy(self) -> Optional[float]:
        if self._occupancy or self._occupancy == 0:
            return self._occupancy
        return None

    @property
    def temperature_factor(self) -> Optional[float]:
        if self._temperature_factor or self._temperature_factor == 0:
            return self._temperature_factor
        return None

    @property
    def _reconstructed_line(self):
        if self.is_het:
            lstart = "HETATM"
        else:
            lstart = "ATOM  "

        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        if len(self.atom_label) == 1:
            atom_label = f" {self.atom_label}  "
        elif len(self.atom_label) == 2:
            if self.atom_label.upper() == self.element.upper():
                atom_label = f"{self.atom_label}  "
            else:
                atom_label = f" {self.atom_label} "
        elif len(self.atom_label) == 3:
            # Undocumented feature of PDBe:
            # If atom name is composite of element and int, padding differs
            ele_len = len(self.element)
            if self.atom_label[:ele_len] == self.element and self.atom_label[ele_len:].isdigit():
                atom_label = "{:>2}{:<2}".format(self.atom_label[:ele_len], self.atom_label[ele_len:])
            else:
                atom_label = "{:>4}".format(self.atom_label)
        elif len(self.atom_label) == 4:
            atom_label = "{:>4}".format(self.atom_label)
        else:
            raise IndexError(f"{self.atom_label} has more than 4 letters!")

        if self._occupancy or self._occupancy == 0:
            occ = "{:0.2f}".format(self._occupancy)
        else:
            occ = " "
        if self._temperature_factor or self._temperature_factor == 0:
            temperature_factor = "{:0.2f}".format(self._temperature_factor)
        else:
            temperature_factor = " "
        fill_values = [lstart,
                       str(self.atom_id),
                       atom_label,
                       self.alt_pos,
                       self.res_name,
                       self.chain,
                       self.res_id,
                       self.icode,
                       "{:0.3f}".format(self.coordinates[0]),
                       "{:0.3f}".format(self.coordinates[1]),
                       "{:0.3f}".format(self.coordinates[2]),
                       occ,
                       temperature_factor,
                       self.segment_id,
                       self.element,
                       self.charge]
        line = "{}{:>5} {:<4}{:>1}{:>3}{:>2}{:>4}{}   {:>8}{:>8}{:>8}{:>6}{:>6}      {:<4}{:>2}{:>2}".format(*fill_values)
        return line

    @classmethod
    def from_pdb_line(cls, line: str, check_line: bool = True):
        """Creates an PDBAtom object from a line in a PDBFile.

        Parameters
        __________
        line: str
            line of PDB file
        check_line: bool
            checks if reconstructed line of atoms is equal to input line. Will raise AssertionError otherwise.

        References
        __________
        www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        """

        if line[:4] == "ATOM":
            is_het = False
        elif line[:6] == "HETATM":
            is_het = True
        else:
            raise ValueError(f"Unknown format: {line[:7]}")
        atom_id = int(line[6:11].lstrip(" "))
        atom_name = line[12:16].strip(" ")
        alt_pos = line[16]
        res_name = line[17:20].strip(" ")
        chain = line[21].strip(" ")
        res_id = int(line[22:26].strip(" "))
        icode = line[26]
        x = float(line[30:38].strip(" "))
        y = float(line[38:46].strip(" "))
        z = float(line[46:54].strip(" "))
        occ = line[54:60].strip(" ")
        if occ:
            occ = float(occ)
        temperature_factor = line[60:66].strip(" ")
        if temperature_factor:
            temperature_factor = float(temperature_factor)
        segment_id = line[72:76].strip(" ")
        element_symbol = line[76:78].strip(" ")
        charge = line[78:80].strip(" ")
        atom = cls(is_het=is_het,
                   atom_id=atom_id,
                   atom_label=atom_name,
                   alt_pos=alt_pos,
                   res_name=res_name,
                   chain=chain,
                   res_id=res_id,
                   icode=icode,
                   coordinates=[x, y, z],
                   occupancy=occ,
                   temperature_factor=temperature_factor,
                   segment_id=segment_id,
                   element=element_symbol,
                   charge=charge)
        if atom.line != line and check_line:
            if atom.line[:len(line)] != line:
                print(atom._reconstructed_line)
                print(line)
                raise AssertionError
        return atom


class PDBResidue:
    METAL_IDS = {"NA", "K", "MG", "CA", "MN", "FE", "FE2", "CO", "NI", "CU", "ZN"}

    def __init__(self,
                 name: str,
                 res_id: int,
                 chain: str,
                 is_hetero: Optional[bool] = None,
                 ):
        """ Object representing a PDB residue

        :param name:
        :param res_id:
        :param chain:
        :param is_hetero:
        """
        self.name: str = name
        self.res_id: int = res_id
        self.chain: str = chain
        self._atoms: List[PDBAtom] = []
        self._het: bool = is_hetero

    def __str__(self):
        return "\t".join([self.chain, self.name, str(self.res_id)])

    def add_atom(self, atom: PDBAtom):
        if atom.chain != self.chain:
            raise ValueError("Atom and Residue must be in the same chain! {}, {}".format(atom.chain, self.chain))
        if atom.res_id != self.res_id:
            raise ValueError("Atom and Residue must be in the same residue-ID! {}, {}".format(atom.res_id, self.res_id))
        self._atoms.append(atom)

    def coordinates(self, no_hydrogen=True) -> List[List[float]]:
        if no_hydrogen:
            return [atom.coordinates for atom in self._atoms if atom.element != "H"]
        else:
            return [atom.coordinates for atom in self._atoms]

    @property
    def atoms(self):
        return self._atoms

    @property
    def block_str(self) -> str:
        return "\n".join([a.line for a in self._atoms])

    @property
    def is_metal(self) -> bool:
        return self.name in PDBResidue.METAL_IDS

    @property
    def is_hetero(self):
        if self._het is None:
            raise ValueError
        else:
            return self._het

    @property
    def text(self) -> str:
        text = []
        for atom in self._atoms:
            text.append(atom.line)
        return "\n".join(text)

    def remove_alternate_positions(self, keep="A"):
        # Checking if alternate positions are available
        for atom in self.atoms:
            if atom.alt_pos != " ":
                break
        else:  # if no breaks
            return self

        # Determininig duplicates
        atom_counter = defaultdict(lambda: 0)
        for atom in self.atoms:
            atom_counter[atom.atom_label] += 1
        atom_duplicates = [atom_label for atom_label, count in atom_counter.items() if count > 1]

        new_atomlist = []
        for atom in self.atoms:
            if atom.atom_label in atom_duplicates:
                if atom.alt_pos == keep:
                    atom.alt_pos = " "
                else:
                    continue
            elif atom.alt_pos != " ":
                warnings.warn(f"Detected alternate position flag without corresponding alt pos for res {self.res_id}")
                atom.alt_pos = " "
            new_atomlist.append(atom)
        self._atoms = new_atomlist
        # Checking completeness
        res_atom_labels = {atom.atom_label for atom in self._atoms}
        lost_atoms = set(atom_counter.keys()) - res_atom_labels
        if lost_atoms:
            print(self.res_id)
            print(atom_counter.keys())
            raise KeyError(f"Lost following atoms: {lost_atoms}")
        return self


class PDBStructure:
    def __init__(self, pdb_string, check_lines=True):
        self.check_lines = check_lines
        self._pdb_string = pdb_string
        self._extract_residues()
        self._extract_links()
        self._extact_connections()

    def _extract_residues(self):
        self._residues: Dict[Tuple[str, int], PDBResidue] = dict()
        self._atoms: Dict[int, PDBAtom] = dict()
        for line in self._pdb_string.split("\n"):
            if line[:4] != "ATOM" and line[:6] != "HETATM":
                continue
            res_name = line[17:20].strip(" ")
            chain = line[20:22].strip(" ")
            res_id = int(line[22:26].strip(" "))
            het = line[:6] == "HETATM"
            if (chain, res_id) not in self._residues:
                self._residues[(chain, res_id)] = PDBResidue(res_name, res_id, chain, het)
            atom = PDBAtom.from_pdb_line(line, self.check_lines)
            self._residues[(chain, res_id)].add_atom(atom)

    def _extract_links(self):
        links = []
        for line in self._pdb_string.split("\n"):
            # TODO Make lines reproducible from properties
            if "LINK   " not in line[:8]:
                continue
            # res = (chain, residue_id)
            res1 = (line[20:22].strip(" "), int(line[22:26].strip(" ")))
            res2 = (line[50:52].strip(" "), int(line[52:56].strip(" ")))
            res1_name = line[16:20].strip(" ")
            res2_name = line[46:50].strip(" ")
            if res1_name == "HOH" or res2_name == "HOH":
                continue
            links.append((res1, res2))
        link_dict = defaultdict(set)
        for r1, r2 in links:
            if r1 not in self._residues:
                raise MissingResidueError(r1)
            if r2 not in self._residues:
                raise MissingResidueError(r2)
            r1_obj = self._residues[r1]
            r2_obj = self._residues[r2]
            link_dict[r1_obj].add(r2_obj)
            link_dict[r2_obj].add(r1_obj)
        self.link_dict: Dict[PDBResidue, Set[PDBResidue]] = dict(link_dict)

    def _extact_connections(self):
        # TODO Make lines reproducible from properties
        atom_dict = self.atom_dict
        inter_res_bonds = defaultdict(set)
        for line in self._pdb_string.split("\n"):
            if line[:6] != "CONECT":
                continue
            line = "{:<80}".format(line)
            atom = int(line[6:11].strip())
            neighbor1 = line[11:16].strip()
            neighbor2 = line[16:21].strip()
            neighbor3 = line[21:26].strip()
            neighbor4 = line[26:31].strip()
            neighbor_list = [neighbor1, neighbor2, neighbor3, neighbor4]
            neighbor_list = [int(neighbor) for neighbor in neighbor_list if neighbor]
            atom_obj = atom_dict[atom]
            res_of_atom = self._residues[(atom_obj.chain, atom_obj.res_id)]
            for neighbor in neighbor_list:
                neighbor_obj = atom_dict[neighbor]
                res_of_neighbor = self._residues[(neighbor_obj.chain, neighbor_obj.res_id)]
                if res_of_atom is res_of_neighbor:
                    inter_res_bonds[atom_obj].add(neighbor_obj)
                    inter_res_bonds[neighbor_obj].add(atom_obj)
                else:
                    if res_of_atom not in self.link_dict:
                        self.link_dict[res_of_atom] = set()
                    self.link_dict[res_of_atom].add(res_of_neighbor)
                    if res_of_neighbor not in self.link_dict:
                        self.link_dict[res_of_neighbor] = set()
                    self.link_dict[res_of_neighbor].add(res_of_atom)
            self._inter_res_bonds: Dict[PDBAtom: Set[PDBAtom]] = {k: v for k, v in inter_res_bonds.items()}

    @property
    def atom_dict(self) -> Dict[int, PDBAtom]:
        atom_dict = dict()
        for res in self._residues.values():
            for atom in res.atoms:
                atom_dict[atom.atom_id] = atom
        return atom_dict

    @property
    def het_res(self) -> List[PDBResidue]:
        return [res for res in self._residues.values() if res.is_hetero]

    @property
    def residues(self) -> list[PDBResidue]:
        return list(self._residues.values())

    @property
    def residue_dict(self) -> Dict[Tuple[str, int], PDBResidue]:
        return self._residues.copy()

    def remove_alternate_positions(self, keep="A"):
        for name, res in self._residues.items():
            res.remove_alternate_positions(keep)
        return self


class PDBFile:
    def __init__(self, pdb_id: str, pdb_content: str, check_lines=True):
        self.pdb_id: str = pdb_id
        self._pdb_file: str = pdb_content
        self.models: List[PDBStructure] = []
        self.check_lines = check_lines
        n_warns = 0
        if self.has_model():
            for model_nr in range(self.max_model()):
                model_nr += 1
                try:
                    self.models.append(self._extract_model(model_nr))
                except MissingResidueError:
                    n_warns += 1
        else:
            self.models.append(PDBStructure(self.pdb_file, self.check_lines))
        if n_warns > 0:
            warnings.warn(f"Error in processing pdb id {self.pdb_id}. {n_warns} models were not load.")
        if not self.models:
            raise ValueError("No valid structure was extracted.")

    @property
    def pdb_file(self) -> str:
        return self._pdb_file

    def __str__(self):
        return self.pdb_id

    def _line_generator(self) -> Generator[str, None, None]:
        for line in self._pdb_file.split("\n"):
            yield line

    def has_model(self) -> bool:
        if "MODEL        1" in self.pdb_file:
            return True
        return False

    def max_model(self):
        model_numbers = []
        for line in self.pdb_file.split("\n"):
            if line[:5] == "MODEL":
                number = line[5:].strip(" ")
                number = int(number)
                model_numbers.append(number)

        return max(model_numbers)

    def _extract_model(self, n=1) -> PDBStructure:
        if self.has_model() is False:
            raise ValueError("This PDB-file does not contain multiple models")
        else:
            out = []
            skip = False
            model_start_line = "MODEL{:>9}".format(n)
            for line in self._line_generator():
                if line[:5] == "MODEL":
                    if line[:14] != model_start_line:
                        skip = True
                if not skip:
                    out.append(line)
                if "ENDMDL" in line:
                    skip = False
            return PDBStructure("\n".join(out))

    def save_to(self, path):
        with open(os.path.join(path, f"{str(self)}.pdb"), "w") as outfile:
            outfile.write(self.pdb_file)

    @classmethod
    def from_file(cls, file_path: str, pdb_id: str, check_lines: bool = True):
        """ Opens file in pdb format and returns a `PDBFile-object`
        :param pdb_id: str
            name or id of file
        :param file_path: str
            path to file
        :param check_lines: bool
            checking if reconstructed atom line equals input line.
            (Some PDF files dont match the conventions, mild errors could be fixed.)
        :return: PDBFile
        """
        with open(file_path) as infile:
            pdb_file = "".join(infile.readlines())

        return cls(pdb_id, pdb_file, check_lines)

    @classmethod
    def from_online(cls, pdb_id: str):
        pdb_id = pdb_id.lower()
        url = f"https://www.ebi.ac.uk/pdbe/entry-files/pdb{pdb_id}.ent"
        return cls(pdb_id, query_website(url).decode('ascii'))

    @classmethod
    def from_file_or_online(cls, pdb_id, folder_path):
        if os.path.isfile(os.path.join(folder_path, f"{pdb_id}.pdb")):
            return cls.from_file(pdb_id, folder_path)
        else:
            return cls.from_online(pdb_id)


def pdb_file_to_dict(pdb_file: PDBFile) -> List[Dict[str, Any]]:
    structure = pdb_file.models[0]
    dict_list = []
    for het_res in structure.het_res:
        if het_res.name == "HOH":
            continue
        ligand_dict = {"pdb_id": pdb_file.pdb_id,
                       "chain": het_res.chain,
                       "res_id": het_res.res_id,
                       "res_name": het_res.name,
                       "connected2protein": False,
                       "connected2ligand": False,
                       "connected2metal": False}

        if het_res in structure.link_dict:
            for connected_res in structure.link_dict[het_res]:
                if connected_res.name == "HOH":
                    continue
                if connected_res.is_metal:
                    ligand_dict["connected2metal"] = True
                elif connected_res.is_hetero:
                    ligand_dict["connected2ligand"] = True
                else:
                    ligand_dict["connected2protein"] = True
        dict_list.append(ligand_dict)

    if not dict_list:
        dict_list = [{"pdb_id": pdb_file.pdb_id,
                      "chain": None,
                      "res_id": None,
                      "res_name": None,
                      "connected2protein": None,
                      "connected2ligand": None,
                      "connected2metal": None,
                      }]
    return dict_list
