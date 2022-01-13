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


def query_website(request, max_trials=10):
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
    def __init__(self,
                 atom_id: int,
                 element: str,
                 res_id: int,
                 chain: str,
                 coordinates: Tuple[float, float, float]):
        self.atom_id = atom_id
        self.symbol = element
        self.res_id = res_id
        self.chain = chain
        self.coordinates = coordinates
        self._line = None

    @property
    def line(self):
        if self._line:
            return self._line
        else:
            raise NotImplementedError("writing text from artificial atoms (not generated from PDB file) is not "
                                      "supported yet!")

    @classmethod
    def from_pdb_line(cls, line):
        atom_id = int(line[6:11].lstrip(" "))
        atom_name = line[11:17].strip(" ")
        res_name = line[17:20].strip(" ")
        chain = line[20:22].strip(" ")
        res_id = int(line[22:26].strip(" "))
        x = float(line[30:38].strip(" "))
        y = float(line[38:46].strip(" "))
        z = float(line[46:54].strip(" "))
        element_symbol = line[76:78].strip(" ")
        atom = cls(atom_id, element_symbol, res_id, chain, (x, y, z))
        atom._line = line
        return atom


class PDBResidue:
    METAL_IDS = {"NA", "K", "MG", "CA", "MN", "FE", "FE2", "CO", "NI", "CU", "ZN"}

    def __init__(self, name: str, res_id: int, chain: str, is_hetero: Optional[bool] = None):
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

    def coordinates(self, no_hydrogen=True) -> List[Tuple[float, float, float]]:
        if no_hydrogen:
            return [atom.coordinates for atom in self._atoms if atom.symbol != "H"]
        else:
            return [atom.coordinates for atom in self._atoms]

    @property
    def atoms(self):
        return self._atoms

    @property
    def block_str(self) -> str:
        return "\n".join([a._line for a in self._atoms])

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


class PDBStructure:
    def __init__(self, pdb_string):
        self._pdb_string = pdb_string
        self._extract_residues()
        self._extract_links()

    def _extract_residues(self):
        self._residues: Dict[Tuple[str, int], PDBResidue] = dict()
        for line in self._pdb_string.split("\n"):
            if line[:4] != "ATOM" and line[:6] != "HETATM":
                continue
            res_name = line[17:20].strip(" ")
            chain = line[20:22].strip(" ")
            res_id = int(line[22:26].strip(" "))
            het = line[:6] == "HETATM"
            if res_name == "HOH":
                het = False
            if (chain, res_id) not in self._residues:
                self._residues[(chain, res_id)] = PDBResidue(res_name, res_id, chain, het)
            self._residues[(chain, res_id)].add_atom(PDBAtom.from_pdb_line(line))

    def _extract_links(self):
        links = []
        for line in self._pdb_string.split("\n"):
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

    @property
    def het_res(self) -> List[PDBResidue]:
        return [res for res in self._residues.values() if res.is_hetero]

    @property
    def residues(self) -> list[PDBResidue]:
        return list(self._residues.values())

    @property
    def residue_dict(self) -> Dict[Tuple[str, int], PDBResidue]:
        return self._residues.copy()


class PDBFile:
    def __init__(self, pdb_id: str, pdb_content: str):
        self.pdb_id: str = pdb_id
        self._pdb_file: str = pdb_content
        self.models: List[PDBStructure] = []
        n_warns = 0
        if self.has_model():
            for model_nr in range(self.max_model()):
                model_nr += 1
                try:
                    self.models.append(self._extract_models(model_nr))
                except MissingResidueError:
                    n_warns += 1
        else:
            self.models.append(PDBStructure(self.pdb_file))
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

    def _extract_models(self, n=1):
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
    def from_file(cls, pdb_id, folder_path):
        with open(os.path.join(folder_path, f"{pdb_id}.pdb")) as infile:
            pdb_file = "".join(infile.readlines())
        return cls(pdb_id, pdb_file)

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
        ligand_dict = {"pdb_id": pdb_file.pdb_id,
                       "chain": het_res.chain,
                       "res_id": het_res.res_id,
                       "res_name": het_res.name,
                       "connected2protein": False,
                       "connected2ligand": False,
                       "connected2metal": False}

        if het_res in structure.link_dict:
            for connected_res in structure.link_dict[het_res]:
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
