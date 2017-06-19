from reformat.reformat import Reformat

if __name__ == '__main__':
    reformat = Reformat()
    reformat.read_input_files()
    reformat.remove_decoy_hits()
    reformat.replace_sample_intensity_column()
    reformat.reformat_peptide_column()
    reformat.reformat_protein_column()
    reformat.write_output_file()
