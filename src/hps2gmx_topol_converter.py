#!/usr/bin/env python3
# hps2gmx_topol_converter written by Xianshi Liu @ Fudan University, 2022. 

import sys
import copy


class Hps2GmxTopol:
    def __init__(self, args:list):
        self.fnhps_file:str = 'hps_topol.txt'
        self.fntopol_file:str = 'topol.top'
        self.parse_args(args)
        self.convert()


    def parse_args(self, args:list):
        arg_count = 1
        while arg_count < len(args):
            if args[arg_count] == '-h':
                self.show_manual()
            elif args[arg_count] == '-hps':
                self.fnhps_file = args[arg_count + 1]
                arg_count += 1
            elif args[arg_count] == '-topol':
                self.fntopol_file = args[arg_count + 1]
                arg_count += 1
            else:
                print("Unknown command line option: " + sys.argv[arg_count])
                quit()

            arg_count += 1


    def convert(self):
        type_list = list()
        type_set = list()
        atom_types = list()
        atom_masses = list()
        atom_charges = list()
        atom_coordinates = list()
        bond_list = list()
        pair_list = list()
        has_charge = False
        next_expect_bond_index = 0
        protein_chain_ter = False
        protein_chain_length = 0
        protein_chain_number = 0
        RNA_chain_ter = False
        RNA_chain_length = 0
        RNA_chain_number = 0

        pbond_length = 0.38
        rbond_length = 0.5
        pbond_number = 0
        rbond_number = 0
        bond_number = 0

        try:
            with open(self.fnhps_file, 'r') as hps_file:
                lines = [line for line in hps_file.readlines() if len(line.strip()) > 0 and line[0] != '#']
                line_index = 0

                while line_index < len(lines):
                    line = lines[line_index]
        
                    if line.strip().split()[0] == "types":
                        types_number = int(line.strip().split()[1])
                        line_index += 1
                        type_list = copy.deepcopy(lines[line_index].strip().split())
                        line_index += 1
        
                    elif line.strip().split()[0] == "atom":
                        atom_number = int(line.strip().split()[1])
                        for atom_index in range(atom_number):
                            line_index += 1
                            atom_types.append(lines[line_index].strip().split()[1])
                            atom_masses.append(float(lines[line_index].strip().split()[2]))
                            atom_charges.append(float(lines[line_index].strip().split()[3]))
                        line_index += 1
        
                    elif line.strip().split()[0] == "pbonds":
                        pbond_number = int(line.strip().split()[1])
                        line_index += 1
        
                    elif line.strip().split()[0] == "rbonds":
                        rbond_number = int(line.strip().split()[1])
                        line_index += 1
        
                    elif line.strip().split()[0] == "bonds":
                        bond_number = int(line.strip().split()[1])
                        for bond_index in range(bond_number):
                            line_index += 1
                            bond_list.append([int(lines[line_index].strip().split()[0]),
                                              int(lines[line_index].strip().split()[1])])
                            if (bond_index < pbond_number 
                                    and protein_chain_ter == False 
                                    and (next_expect_bond_index != bond_list[-1][0] or bond_index == pbond_number - 1)):
                                if (bond_index == pbond_number - 1):
                                    protein_chain_length = bond_list[-1][0] + 2
                                else:
                                    protein_chain_length = bond_list[-1][0]
                                protein_chain_number = int(pbond_number / (protein_chain_length - 1))
                                #print(protein_chain_length)
                                #print(protein_chain_number)
                                protein_chain_ter = True
                            if (bond_index > pbond_number 
                                    and RNA_chain_ter == False 
                                    and (next_expect_bond_index != bond_list[-1][0] or bond_index == bond_number - 1)):
                                if (bond_index == bond_number - 1):
                                    RNA_chain_length = int(bond_list[-1][0] - protein_chain_length * protein_chain_number) + 2
                                else:
                                    RNA_chain_length = int(bond_list[-1][0] - protein_chain_length * protein_chain_number)
                                RNA_chain_number = int(rbond_number / (RNA_chain_length - 1))
                                RNA_chain_ter = True
                            next_expect_bond_index = bond_list[-1][1]
                        line_index += 1
        
                    else:
                        print("ERROR: Cannot process input parameter file line:")
                        print(lines[line_index])
                        quit()

        except FileNotFoundError:
            print(f"ERROR: Cannot load paramter file: {self.fnhps_file}")
            quit()


        try:
            with open(self.fntopol_file, 'w') as topol_file:
                topol_file.write("; HPS from xsliu\n\n")
                topol_file.write("[ defaults ]\n")
                topol_file.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
                topol_file.write("1               2               yes              0.5     0.8333\n\n")
                topol_file.write("[ atomtypes ]\n")
                topol_file.write(";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n")

                for atom_type in type_list:
                    if atom_type in aminoacids:
                        if aminoacids[atom_type] not in type_set:
                            topol_file.write(f"\t{aminoacids[atom_type]}\t{aminoacids[atom_type]}\t0\t0\tA\t1\t1\n")
                            type_set.append(aminoacids[atom_type])
                    elif atom_type in nucleicacids:
                        if nucleicacids[atom_type] not in type_set:
                            topol_file.write(f"\t{nucleicacids[atom_type]}\t{nucleicacids[atom_type]}\t0\t0\tA\t1\t1\n")
                            type_set.append(nucleicacids[atom_type])
                    else:
                        print(f"ERROR: Cannot parse the type: {atom_type}")
                        quit()


                if pbond_number > 0:
                    topol_file.write("\n[ moleculetype ]\n")
                    topol_file.write(";name            nrexcl\n")
                    topol_file.write(" protein            1\n\n")
                    topol_file.write("[ atoms ]\n")
                    topol_file.write(";   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type\n")
    
                    for atom_index in range(protein_chain_length):
                        if atom_types[atom_index] in aminoacids:
                            topol_file.write(f"{atom_index + 1}\t{aminoacids[atom_types[atom_index]]}\t{atom_index + 1}" +
                                f"\t{aminoacids[atom_types[atom_index]]}\t{aminoacids[atom_types[atom_index]]}\t{atom_index + 1}" +
                                f"\t{atom_charges[atom_index]}\t{atom_masses[atom_index]}\n")
                        else:
                            print(f"ERROR: Cannot parse the atom: {atom_types[atom_index]}")
                            quit()
    
    
    
                    topol_file.write("\n[ bonds ]\n")
                    topol_file.write(";   ai     aj funct   ARG             k\n")
    
                    for bond_index in range(protein_chain_length - 1):
                        topol_file.write("\t".join(str(i + 1) for i in bond_list[bond_index]))
                        topol_file.write(f"\t1\t{pbond_length}\t1\n")


                if rbond_number > 0:
                    topol_file.write("\n[ moleculetype ]\n")
                    topol_file.write(";name            nrexcl\n")
                    topol_file.write(" RNA                1\n\n")
                    topol_file.write("[ atoms ]\n")
                    topol_file.write(";   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type\n")
    
                    for atom_index_r in range(RNA_chain_length):
                        atom_index = int(atom_index_r + protein_chain_length * protein_chain_number)
                        if atom_types[atom_index] in nucleicacids:
                            topol_file.write(f"{atom_index_r + 1}\t{nucleicacids[atom_types[atom_index]]}\t{atom_index_r + 1}" +
                                f"\t{nucleicacids[atom_types[atom_index]]}\t{nucleicacids[atom_types[atom_index]]}\t{atom_index_r + 1}" +
                                f"\t{atom_charges[atom_index]}\t{atom_masses[atom_index]}\n")
                        else:
                            print(f"ERROR: Cannot parse the atom: {atom_types[atom_index]}")
                            quit()
    
    
    
                    topol_file.write("\n[ bonds ]\n")
                    topol_file.write(";   ai     aj funct   ARG             k\n")
    
                    for bond_index_r in range(RNA_chain_length - 1):
                        bond_index = bond_index_r + pbond_number
                        topol_file.write("\t".join(str(i + 1 - protein_chain_length * protein_chain_number) for i in bond_list[bond_index]))
                        topol_file.write(f"\t1\t{rbond_length}\t1\n")



                topol_file.write("\n\n[ system ]\n")
                topol_file.write("HPS\n\n")
                topol_file.write("[ molecules ]\n")
                topol_file.write("; Compound        nmols\n")
                if pbond_number > 0:
                    topol_file.write(f" protein            {protein_chain_number} \n")
                if rbond_number > 0:
                    topol_file.write(f" RNA                {RNA_chain_number} \n")

            print("Successfully converted!!!")

        except FileNotFoundError:
            print(f"ERROR: Cannot open file: {self.fntopol_file}")
            quit()



    def show_manual(self):
        manual = (
            "Program: hps2gmx_topol_converter written by Xianshi Liu @ Fudan\n" +
            "This program convert a HOOMD-Blue input file to Gromacs topology file.\n" +
            "Parameters for file options:\n" +
            "      -hps   [hps_topol.txt]   : Topology information file\n" +
            "      -topol [topol.top]       : GROMACS output topology file\n" +
            "      -h                       : Show this manual\n")

        print(manual)
        quit()
        

if __name__ == '__main__':
    Hps2GmxTopol(sys.argv)


