# PDBPython
The scope of this module is downloading or reading PDB-Files and subsequently parsing them.

## Disclaimer
This project is only intended for my personal use and may change significantly over time. You are free to use this code
at own risk.

## Classes
### PDBFile
The class PDBFile is initialized with `pdb_id` (str, an identifier for the file) and `pdb_content` (str, the content of 
a PDB-file). In default `check_lines` (bool) is set to `True`, which will raise an `AssertionError` if the reconstructed
line for an atom differs from the original line. In some cases PDB files created by third party software may differ
slightly in padding, which will still raise an error. After reviewing differences you can suppress the error  with
setting `check_lines` to `False`.
#### Constructing a PDBFile object from online database or file.
The class PDBFile can also be constructed from a file:
```python
pdbfile_obj = PDBFile.from_file("./path/to/file.pdb", "PDB_ID_of_File")
```
Or directly from [PDBe](https://www.ebi.ac.uk/pdbe/)
```python
pdbfile_obj = PDBFile.from_online("PDB_ID")
```
### PDBStructure
Most PDB files contain only one protein model, however there are PDB files which contain multiple protein conformations
(mainly NMR-resolved structures).
Models are herein called `PDBStructure` and can be accessed via the PDBFIle object:
```python
pdbfile_obj = PDBFile.from_online("PDB_ID")
model_0 = pdbfile_obj.model[0]
```
On this note: My naming conventions with model and structure may be not very consistent. This is a #TODO for future
refactorings.  

### PDBResidue
Each PDBStructure consists of residudes represented by `PDBResidue`-objects. Which themselves consist of `PDBAtom`-
objects.  
Residues can be accessed via:
```python
pdbfile_obj = PDBFile.from_online("PDB_ID")
model_0 = pdbfile_obj.model[0]
residue = model_0.residues
```
Residues specified by chain and residue-ID can also be accessed via a dictionary:
```python
chain = "A"
res_id = 1
residue_A_1 = model_0.residue_dict[(chain, res_id)]
```
For some atoms in a residue alternative positions are given. When specific atom positions are required other locations 
can be removed via:
```python
residue.remove_alternate_positions(keep="A")
```
`keep` specifies the alternative position to keep. The default value is `"A"`.  
To apply this to all residues of a PDBStructure, call:
```python
model_0.remove_alternate_positions(keep="A")
```