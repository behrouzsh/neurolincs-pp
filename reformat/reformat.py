import pandas as pd
class Reformat(object):
    """reformat class is to read the input file and other auxiliary files and reformat it in the
     desired format.
    """


    def __init__(self):
        self.motif_and_modification_list = []
        self.prosite_result = []
        self.uniprot_result = []
        self.psimod_result = []
        self.unique_protein = []
        self.unique_genes = []
        self.network = []
        self.pathway_list = []
        self.gene_to_pathway_list = []

    def read_input_files(self):
        with open('resource/Protein_iPSC_Metdata_Master_modified.xlsx', "r") as xls_data_file:
            self._metadata = pd.read_excel(xls_data_file)

        with open('resource/Renaming_neuroLINCS.xlsx', "r") as xls_data_file:
            self._metadata2 = pd.read_excel(xls_data_file, header = None)

        # This file is coming from www.nature.com/
        # Just search for srep07896-s9.xls in google

        with open('resource/srep07896-s9.xls', "r") as xls_data_file:
            self._uniModdata = pd.read_excel(xls_data_file)


        xls = pd.ExcelFile('input/outputmatrix_OpenSWATH.xlsx')
        self._input_data = xls.parse('IPSC_DDA_6600_QE_combined_Canon')

        return

    def remove_decoy_hits(self):

        self._input_data = self._input_data[self._input_data['Peptide'].str.contains("DECOY") == False]

        return


    def replace_sample_intensity_column(self):
        # Get the header names
        header_names = list(self._input_data)

        renamed_columns = self.find_new_sample_names2(header_names)

        self._input_data.rename(columns = renamed_columns, inplace = True)

        return

    def find_new_sample_names(self, pre_sample_names):
        """This function Replace sample 'intensity_' column headers with official sample generation identifiers
            @param: pre_sample_names The column names of the input data
            """
        new_sample_names = {}
        Intensity_MS_File_Name = self._metadata['Intensity_MS_File_Name']

        IP_Center_Specific_ID = self._metadata['IP_Center_Specific_ID']
        for sample in range (0, len(Intensity_MS_File_Name)):
            for pre_sample in pre_sample_names:

                stripped = str(Intensity_MS_File_Name[sample])[5:].strip()
                if stripped in str(pre_sample):
                    new_sample_names[str(pre_sample)] = str(IP_Center_Specific_ID[sample])

        return new_sample_names

    def find_new_sample_names2(self, pre_sample_names):
        """This function Replace sample 'intensity_' column headers with official sample generation identifiers
            @param: pre_sample_names The column names of the input data
            """
        new_sample_names = {}

        Intensity_MS_Sample_Name = self._metadata2.loc[0,:]


        IP_Center_Specific_ID = self._metadata2.loc[1,:]

        for sample in range (0, len(Intensity_MS_Sample_Name)):
            for pre_sample in pre_sample_names:
                if str(Intensity_MS_Sample_Name[sample]) == str(pre_sample):
                    new_sample_names[str(pre_sample)] = str(IP_Center_Specific_ID[sample])

        return new_sample_names

    def reformat_peptide_column(self):
        unimod_peptide = []
        amu_peptide = []
        peptide_sequence = []
        charge = []
        peptide_col = self._input_data['Peptide']
        for peptide in peptide_col:
            splitted_by_underscore = str(peptide).split("_")
            pep = splitted_by_underscore[1]
            unimod_peptide.append(pep)
            charge.append(splitted_by_underscore[2])
            if "(" in pep:
                unimod = pep[pep.find("(") + 1:pep.find(")")]
                sequence = pep[0:pep.find("(")] + pep[pep.find(")") + 1:]
            else:
                unimod = ""
                sequence = pep

            peptide_sequence.append(sequence)

            # Find amu increase
            if unimod != "":
                amu_increase = self.find_amu_increase(unimod)
                peptide_plus_amu = pep[0:pep.find("(")] + "(" + str(amu_increase) + ")" + pep[pep.find(")") + 1:]
            else:
                peptide_plus_amu = pep
            amu_peptide.append(peptide_plus_amu)



        self._input_data['Unimod_peptide'] = unimod_peptide
        self._input_data['Amu_peptide'] = amu_peptide
        self._input_data['Peptide_sequence'] = peptide_sequence
        self._input_data['Charge(z)'] = charge
        return

    def find_amu_increase(self, unimod):
        """This function finds the amu increase of the unimod
            @param: unimod The unimod of PTM
            """

        unimod_access = int(unimod.split(":")[1])

        unimod_loc = self._uniModdata[self._uniModdata['UniMod: Accession #'] == unimod_access].index.tolist()


        amu_mass = float(self._uniModdata['UniMod: Monoisotopic Mass (not from UniMod when red)'][unimod_loc[0]])


        return amu_mass


    def reformat_protein_column(self):
        protein_num = []
        uniProtKB = []
        gene_organism = []
        protein_col = self._input_data['Protein']
        for protein in protein_col:
            splitted_by_pipe = str(protein).split("|")
            prot_num = splitted_by_pipe[0].split("/")[0]
            protein_num.append(prot_num)

            uniProtKB_str = ""
            gene_organism_str = ""
            for iter in range(1, len(splitted_by_pipe),2):
                uniProtKB_str = uniProtKB_str + " " + splitted_by_pipe[iter]

            for iter in range(2, len(splitted_by_pipe),2):
                gene_organism_str = gene_organism_str + " " + splitted_by_pipe[iter].split("/")[0]

            uniProtKB.append(uniProtKB_str)
            gene_organism.append(gene_organism_str)

        self._input_data['Number of proteins mapped'] = protein_num
        self._input_data['UniProtKB'] = uniProtKB
        self._input_data['Gene_organism'] = gene_organism

        return


    def write_output_file(self):
        self._input_data.to_csv('output/outputmatrix_OpenSWATH.csv', sep=',')





