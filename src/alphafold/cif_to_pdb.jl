using BioStructures

@assert length(ARGS) == 2
in_cif_file = ARGS[1]
out_pdb_file = ARGS[2]

#extension = split(in_cif_file, ".")[end]
#out_pdb_file = in_cif_file[1:end-length(extension)-1] * ".pdb"

mmcif_dict = MMCIFDict(in_cif_file)
mmcif_dict["_atom_site.auth_atom_id"] = mmcif_dict["_atom_site.label_atom_id"]
mmcif_dict["_atom_site.auth_comp_id"] = mmcif_dict["_atom_site.label_comp_id"]
mmcif_dict["_atom_site.pdbx_formal_charge"] = fill("?", length(mmcif_dict["_atom_site.label_atom_id"]))
ms = MolecularStructure(mmcif_dict)
writepdb(out_pdb_file, ms)
