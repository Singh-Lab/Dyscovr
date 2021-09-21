#########################################################################################################
### AGGREGATE MAF FILES ###
### Sara Camilli, July, 2020 ###
#########################################################################################################

import sys
import os
import sys
import gzip
import math

# This file takes in a bunch of MAF files and aggregates them into a single output MAF file
def main():
    output_maf_filename = "TCGA.BRCA.Aggregated.Muse.maf"

    input_maf_paths = ["0cdf3c70-ad58-462d-b6ba-5004b26c618e/TCGA.LAML.muse.0cdf3c70-ad58-462d-b6ba-5004b26c618e.DR-10.0.somatic.maf",
     "15b448f9-56f8-4ba5-8959-674312ee9cfa/TCGA.UCS.muse.15b448f9-56f8-4ba5-8959-674312ee9cfa.DR-10.0.somatic.maf",
     "21a77f28-8512-4224-9987-0eef620511e3/TCGA.THYM.muse.21a77f28-8512-4224-9987-0eef620511e3.DR-10.0.somatic.maf",
     "225a3a8f-27fd-4d9f-b8af-4eadbae9304e/TCGA.HNSC.muse.225a3a8f-27fd-4d9f-b8af-4eadbae9304e.DR-10.0.somatic.maf",
     "27330980-b14d-4d4d-bd65-5eb19c597c9c/TCGA.LGG.muse.27330980-b14d-4d4d-bd65-5eb19c597c9c.DR-10.0.somatic.maf",
     "2a332e6f-5a76-43cb-aa78-4cd746891b1b/TCGA.SARC.muse.2a332e6f-5a76-43cb-aa78-4cd746891b1b.DR-10.0.somatic.maf",
     "39f0a058-6ce5-4bbc-bcc8-57d4cd14066f/TCGA.CESC.muse.39f0a058-6ce5-4bbc-bcc8-57d4cd14066f.DR-10.0.somatic.maf",
     "3f2602b5-01dc-49f0-8393-33a4ca2c64f2/TCGA.KICH.muse.3f2602b5-01dc-49f0-8393-33a4ca2c64f2.DR-10.0.somatic.maf",
     "51423d79-e9c5-4c4d-b12c-99c1338dbd43/TCGA.OV.muse.51423d79-e9c5-4c4d-b12c-99c1338dbd43.DR-10.0.somatic.maf",
     "58390670-c78f-4baf-b478-90ee01d15158/TCGA.STAD.muse.58390670-c78f-4baf-b478-90ee01d15158.DR-10.0.somatic.maf",
     "59531694-1281-4053-9516-15e146022fcb/TCGA.PRAD.muse.59531694-1281-4053-9516-15e146022fcb.DR-10.0.somatic.maf",
     "59a84472-27d4-497c-8f37-8bc447ff9374/TCGA.GBM.muse.59a84472-27d4-497c-8f37-8bc447ff9374.DR-10.0.somatic.maf",
     "6d0a808a-bdd7-40b0-94d1-168e73ad5bdb/TCGA.MESO.muse.6d0a808a-bdd7-40b0-94d1-168e73ad5bdb.DR-10.0.somatic.maf",
     "6f5cde97-d259-414f-8122-6d0d66f49b74/TCGA.LUAD.muse.6f5cde97-d259-414f-8122-6d0d66f49b74.DR-10.0.somatic.maf",
     "70cb1255-ec99-4c08-b482-415f8375be3f/TCGA.COAD.muse.70cb1255-ec99-4c08-b482-415f8375be3f.DR-10.0.somatic.maf",
     "716f5244-9e13-4375-8c30-147bf9a02bb4/TCGA.KIRC.muse.716f5244-9e13-4375-8c30-147bf9a02bb4.DR-10.0.somatic.maf",
     "7837f512-39c2-4c7c-a338-43e51a54701c/TCGA.CHOL.muse.7837f512-39c2-4c7c-a338-43e51a54701c.DR-10.0.somatic.maf",
     "7ed5eb26-d52c-45c5-bfa1-e39f963007a8/TCGA.BLCA.muse.7ed5eb26-d52c-45c5-bfa1-e39f963007a8.DR-10.0.somatic.maf",
     "93c525cc-655c-4c1c-b590-18d851473f68/TCGA.PAAD.muse.93c525cc-655c-4c1c-b590-18d851473f68.DR-10.0.somatic.maf",
     "97c1698f-d4a0-4bb0-949b-92e59011f438/TCGA.KIRP.muse.97c1698f-d4a0-4bb0-949b-92e59011f438.DR-10.0.somatic.maf",
     "9d298cc7-71d2-43f5-835d-d7ef7a43bb11/TCGA.ACC.muse.9d298cc7-71d2-43f5-835d-d7ef7a43bb11.DR-10.0.somatic.maf",
     "ad676dc8-bc24-440f-bedd-b0885891fb61/TCGA.LUSC.muse.ad676dc8-bc24-440f-bedd-b0885891fb61.DR-10.0.somatic.maf",
     "b267fee6-e9a0-4d50-88fb-79eb8da1e8dd/TCGA.ESCA.muse.b267fee6-e9a0-4d50-88fb-79eb8da1e8dd.DR-10.0.somatic.maf",
     "b8ca5856-9819-459c-87c5-94e91aca4032/TCGA.BRCA.muse.b8ca5856-9819-459c-87c5-94e91aca4032.DR-10.0.somatic.maf",
     "b9d7ec0e-afac-41ee-bd08-0acbfa0e3858/TCGA.LIHC.muse.b9d7ec0e-afac-41ee-bd08-0acbfa0e3858.DR-10.0.somatic.maf",
     "d12371d7-18ff-4105-a4a0-59de52b82805/TCGA.UCEC.muse.d12371d7-18ff-4105-a4a0-59de52b82805.DR-10.0.somatic.maf",
     "da39b25a-2c4b-413c-a95e-0d65c07da573/TCGA.TGCT.muse.da39b25a-2c4b-413c-a95e-0d65c07da573.DR-10.0.somatic.maf",
     "dcf01652-1420-4189-b087-b968692c6208/TCGA.THCA.muse.dcf01652-1420-4189-b087-b968692c6208.DR-10.0.somatic.maf",
     "e433a47f-7281-42aa-8bba-0affd1eb6df0/TCGA.SKCM.muse.e433a47f-7281-42aa-8bba-0affd1eb6df0.DR-10.0.somatic.maf",
     "ec8ec3ad-f08d-46eb-9571-42806e304b37/TCGA.READ.muse.ec8ec3ad-f08d-46eb-9571-42806e304b37.DR-10.0.somatic.maf",
     "eea22323-9e5b-45e3-90c2-8aa7d07da175/TCGA.DLBC.muse.eea22323-9e5b-45e3-90c2-8aa7d07da175.DR-10.0.somatic.maf",
     "f8dabf59-3ae8-4358-8dd6-5c1df76f7a12/TCGA.PCPG.muse.f8dabf59-3ae8-4358-8dd6-5c1df76f7a12.DR-10.0.somatic.maf"]

    create_aggregate_maf_files(input_maf_paths, output_maf_filename)

    # TODO: Write code to subset the output MAF to only include patients that overlap between all data types

# Function borrowed and modified from Shilpa
def create_aggregate_maf_files(input_mafs, output_maf):
    """
    :param input_mafs: set or list of full paths to maf files to be included in the output maf file
    :param output_maf: full path to an "output" maf file to write aggregate mutations to
    :return: none; this step might require a LOT of memory!
    """

    # (1) process JUST the headers from each of the input mafs to make sure we combine columns correctly
    comments = ['# All mutation data from the following files combined in ' + output_maf + ':']
    master_header = []
    for maf_file in input_mafs:
        if not os.path.isfile(maf_file):
            sys.stderr.write('No maf file ' + maf_file + '\n')
            continue

        current_header = None
        maf_handle = gzip.open(maf_file, 'rt') if maf_file.endswith('gz') else open(maf_file)
        for maf_line in maf_handle:
            if maf_line.startswith('#'):
                continue
            current_header = maf_line[:-1].split('\t')
            break
        maf_handle.close()

        # add all elements of the new header into the master_header
        new_elements = [elem for elem in current_header if elem not in master_header]
        master_header.extend(new_elements)

        # keep track of a running comment header for the output file
        comments.append('# ' + maf_file)

    # (2) NOW we are ready to combine all the maf files
    if len(master_header) < 1:
        sys.stderr.write('No input maf files to be written to ' + output_maf + '!\n')
        return None

    master_handle = gzip.open(output_maf, 'wt') if output_maf.endswith('gz') else open(output_maf, 'w')
    master_handle.write('\n'.join(comments) + '\n')
    master_handle.write('\t'.join(master_header) + '\n')

    seen = set()  # keep track of mutations that we have already seen (to avoid duplicates as we go)
    for maf_file in input_mafs:
        if not os.path.isfile(maf_file):
            continue

        current_header = None
        maf_handle = gzip.open(maf_file, 'rt') if maf_file.endswith('gz') else open(maf_file)
        for maf_line in maf_handle:
            if maf_line.startswith('#'):
                continue
            elif not current_header:
                current_header = maf_line[:-1].split('\t')
                continue
            cvals = maf_line[:-1].split('\t')
            ovals = [cvals[current_header.index(elem)] if elem in current_header else '' for elem in master_header]
            curr_mut = [cvals[current_header.index(elem)] for elem in
                        ['NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Tumor_Sample_Barcode',
                         'Reference_Allele']]
            curr_mut.append(','.join(sorted(list({cvals[current_header.index('Tumor_Seq_Allele1')],
                                                  cvals[current_header.index('Tumor_Seq_Allele2')]}))))

            if tuple(curr_mut) not in seen:  # if this is a unique mutation that has not yet been seen
                master_handle.write('\t'.join(ovals) + '\n')
                seen.add(tuple(curr_mut))
            else:
                sys.stderr.write('Duplicate mutation observed in ' + maf_file + ':\n' + maf_line + '\n')
        maf_handle.close()
    master_handle.close()
    return output_maf

if __name__ == "__main__":
    # execute only if run as a script
    main()