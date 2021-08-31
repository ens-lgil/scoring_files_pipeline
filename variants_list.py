import os
import argparse
import pandas as pd
import gzip


remap_header = {
    'PGS ID': 'pgs_id',
    'PGS Name': 'pgs_name',
    'Reported Trait': 'trait_reported',
    'Original Genome Build': 'genome_build',
    'Number of Variants': 'variants_number',
    'PGP ID': 'pgp_id',
    'Citation': 'citation',
    'LICENSE': 'pgs_license',
    # Harmonization related
    'HmPOS Build': 'HmPOS_build',
    'HmPOS Date':'HmPOS_date',
    'HmVCF Reference': 'HmVCF_ref',
    'HmVCF Date': 'HmVCF_date'
} # ToDo remove once Scoring File headers are fixed


def read_scorefile(loc_scorefile):
    """Loads PGS Catalog Scoring file and parses the header into a dictionary"""
    if loc_scorefile.endswith('.gz'):
        f = gzip.open(loc_scorefile,'rt')
    else:
        f = open(loc_scorefile, 'rt')

    header = {}
    lastline = '#'
    while lastline.startswith('#'):
        lastline = f.readline()
        line = lastline.strip()
        if line.startswith('#'):
            if '=' in line:
                line = line[1:].split('=')
                field, val = [x.strip() for x in line]
                header[remap_header[field]] = val # ToDo change once Scoring File headers are fixed
    f.close()

    if ('genome_build' in header) and (header['genome_build'] == 'NR'):
        header['genome_build'] = None

    df_scoring = pd.read_table(loc_scorefile, float_precision='round_trip', comment='#',
                               dtype = {'chr_position' : 'object',
                                        'chr_name' : 'object',
                                        'hm_chr': 'object',
                                        'hm_pos': 'object'})

    # Make sure certain columns maintain specific datatypes
    if 'reference_allele' in df_scoring.columns:
        df_scoring = df_scoring.rename(columns={"reference_allele": "other_allele"})

    # for poscol in ['chr_position', 'hm_pos']:
    #     if poscol in df_scoring.columns:
    #         df_scoring.loc[df_scoring[poscol].isnull() == False, poscol] = [conv2int(x) for x in df_scoring.loc[df_scoring[poscol].isnull() == False, poscol]]

    return(header, df_scoring)


def clean_rsIDs(raw_rslist):
    """Takes a list of values, removes anything that doesn't look like an rsID and splits any variants that
    are haplotypes, combinations, or interactions"""
    cln_rslist = set()
    for x in raw_rslist:
        if type(x) is str and x.startswith('rs'):
            if '_x_' in x:
                x = [y.strip() for y in x.split('_x_')]
            elif ';' in x:
                x = [y.strip() for y in x.split(';')]
            elif ',' in x:
                x = [y.strip() for y in x.split(',')]
            else:
                cln_rslist.add(x)

            if type(x) == list:
                for i in x:
                    if i.startswith('rs'):
                        cln_rslist.add(i)
    return(list(cln_rslist))


def read_coord_vars_file(coord_file):
    rsIDs_list = set()
    with open(coord_file) as file:
        for line in file:
            line_components = line.split('\t')
            rsID_1 = line_components[0]
            rsID_2 = line_components[1]
            rsIDs_list.add(rsID_1)
            rsIDs_list.add(rsID_2)
    return rsIDs_list


def main():

    argparser = argparse.ArgumentParser()
    argparser.add_argument("--scores_ids", help='List of scores IDs', required=True, metavar='PGS_IDS')
    argparser.add_argument("--scores_dir", help='Directory hosting the scoring files', required=True, metavar='PGS_DIR')
    argparser.add_argument("--var_file", help='Path to the variants output file', required=True, metavar='VAR_FILE')
    argparser.add_argument("--coord_file", help='Path to the variants file with coordinates', required=True, metavar='VAR_COORD_FILE')

    args = argparser.parse_args()

    scores_dir = args.scores_dir
    var_file = args.var_file
    coord_file = args.coord_file

    if not os.path.isdir(scores_dir):
        print("Directory '"+scores_dir+"' can't be found")
        exit(1)

    if not os.path.isfile(coord_file):
        print("File '"+coord_file+"' can't be found")
        exit(1)

    
    existing_rsIDs_with_coords = read_coord_vars_file(coord_file)

    var_list = set()

    for score_id in args.scores_ids.split(','):
        scorefile = scores_dir+'/'+score_id+'.txt.gz'
        if not os.path.isfile(scorefile):
            print("File '"+scorefile+"' can't be found")
            exit(1)
        
        header, df_scoring = read_scorefile(scorefile)
        rsIDs = clean_rsIDs(list(df_scoring['rsID']))

        # print(f'> {score_id}: {len(rsIDs)}')

        for rsID in rsIDs:
            if not rsID in var_list  and not rsID in existing_rsIDs_with_coords:
                var_list.add(rsID)
    
    print(f'Variants: {len(var_list)}')

    # Create "variants" directory if it doesn't exist
    var_dir = os.path.dirname(var_file)
    if not os.path.isdir(var_dir):
        os.makedirs(var_dir)

    # Print the list of variants in a file
    file_out = open(var_file, 'w')
    for var_id in var_list:
        file_out.write(f'{var_id}\n')
    file_out.close()


if __name__ == '__main__':
    main()