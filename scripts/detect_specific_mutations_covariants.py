import pandas as pd
import os
import datetime

def clean_raw_muts(muts):
    """Clean raw mutation strings from freyja output"""
    output = []
    for m in muts.split(" "):
        if ":" not in m or "*" in m:
            continue
        if ":INS" in m:
            insertion_size = len(m.split(")(")[0].split(",'")[1])
            if insertion_size % 3 != 0:
                continue
            output.append(m)
            continue
        if ":DEL" in m:
            deletion_size = m.split(")(")[0].split(",")[1]
            if int(deletion_size) % 3 != 0:
                continue
            output.append(m)
        else:
            output.append(m)
    if len(output) == 0:
        return []
    return output


def parse_aa_muts(muts):
    output = []
    for m in muts:
        if ":DEL" in m:
            output.append(m.split(")(")[1][:-1])
        elif ":INS" in m:
            output.append(m.split(")(")[1][:-1])
        else:
            output.append(m.split("(")[1][:-1])
    return list(set(output))


def parse_covariants(covariants_dir, metadata_file, sample_id):
    """Parse freyja covariants output files, aggregate into one dataframe"""

    agg_covariants = pd.DataFrame()
    for file in os.listdir(covariants_dir):
        df = pd.read_csv(f'{covariants_dir}/{file}', sep="\t")
        try:
            df["Covariants"] = df["Covariants"].apply(clean_raw_muts)
        except:
            continue
        df["aa_muts"] = df["Covariants"].apply(parse_aa_muts)
        df = df[df["aa_muts"].apply(len) > 0]
        df["sample"] = file
        agg_covariants = pd.concat([agg_covariants, df])

    agg_covariants['sample'] = agg_covariants['sample'].apply(lambda x: x.replace('_','-').replace('ENV-','').split('-S')[0])

    # Merge metadata with covariants
    if metadata_file is not None:
        if metadata_file.endswith(".csv"):
            times_df = pd.read_csv(metadata_file)
        else:
            times_df = pd.read_csv(metadata_file, sep="\t")

        times_df['Sequence_ID'] = times_df['Sequence_ID'].apply(lambda x:x.replace('_','-').replace('ENV-','').split('.')[0])
        times_df['LabNumber'] = times_df['LabNumber'].apply(lambda x:x.replace('ENV-',''))
        times_df = times_df.drop_duplicates()
        times_df['SampleCollectionDate'] = pd.to_datetime(times_df['SampleCollectionDate'])
        times_df = times_df.loc[~times_df.duplicated(keep='first')]
        daysIncluded=180
        cut_date = pd.to_datetime(times_df['SampleCollectionDate'].max()-datetime.timedelta(days=daysIncluded))
        times_df = times_df[times_df['SampleCollectionDate']>=cut_date]

        agg_covariants[sample_id] = agg_covariants['sample'].apply(lambda x: x.split(".")[0])
        agg_covariants = pd.merge(agg_covariants, times_df, on=sample_id, how="left")

    return agg_covariants

def main():
    covariants_dir = '../covariants/'
    agg_covariants = parse_covariants(covariants_dir, '../sample_metadata.csv', 'LabNumber')

    # check for mutations of interest.
    query_muts = ['S:T95I','S:I101T','S:N164K','S:S172F','S:DEL136/147']
    agg_covariants = agg_covariants[agg_covariants['aa_muts'].apply(lambda x: any(mut in x for mut in query_muts))]
    agg_covariants.to_csv('covariants_detects.csv', index=False)

if __name__ == "__main__":
    main()