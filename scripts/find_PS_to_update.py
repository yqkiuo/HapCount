import pandas as pd
import argparse
import csv

def process_vcf(shapeit_vcf, whatshap_vcf, output_file):
    # Read input VCF files
    shp = pd.read_csv(shapeit_vcf, sep='\t', comment='#', header=None)
    whp = pd.read_csv(whatshap_vcf, sep='\t', comment='#', header=None)

    # Rename columns
    columns = 'Chromosome POS ID REF ALT Score FILTER INFO FORMAT Sample'.split()
    whp = whp.rename(columns={i: x for i, x in enumerate(columns)})
    shp = shp.rename(columns={i: x for i, x in enumerate(columns)})

    # Process shapeit VCF
    shpt = shp.loc[:, 'Chromosome POS REF ALT Sample'.split()].rename(columns={'Sample': 'ShapeitGT'})

    # Process whatshap VCF
    whp = whp.assign(whatshpGT=[x.split(':')[0] for x in whp.Sample])
    whp = whp.assign(Chromosome=whp.Chromosome.astype(str))
    shpt = shpt.assign(Chromosome=[str(x[3:]) for x in shp.Chromosome])
    whp = whp.assign(PS=[x.split(':')[-1] for x in whp.Sample])

    # Merge data
    shpm = shpt.merge(whp.loc[:, 'Chromosome POS REF ALT whatshpGT PS'.split()], on='Chromosome POS REF ALT'.split())
    shpm = shpm.assign(swapped_by_shapeit=shpm.ShapeitGT == shpm.whatshpGT)

    # Identify strange PS values
    tmp = shpm.groupby('PS').agg({'swapped_by_shapeit': 'nunique'})['swapped_by_shapeit']
    strange_PS = set(tmp[tmp == 2].index)

    if len(list(strange_PS))==0:
            print("No inconsistant PS tags found")
            whp = pd.read_csv(whatshap_vcf, sep='\t', comment='#', header=None)

            # Rename columns
            columns = 'Chromosome POS ID REF ALT Score FILTER INFO FORMAT Sample'.split()
            whp = whp.rename(columns={i: x for i, x in enumerate(columns)})
            whp = whp.assign(Chromosome = whp.Chromosome.astype(str))
            whp.to_csv(output_file, sep='\t', index=False)
            print(f"Output saved to {output_file}")
            return 0






    print("Managing inconsistant PS tags")
    inconst = len(list(strange_PS))
    all_PS = whp.PS.nunique()
    print(f'{inconst} out of {all_PS} PS tags are inconsistant')
    shpm = shpm.assign(PS_copy=shpm.PS)

    # Function to process dataframe
    def process_df(tmp_o):
        tmp = tmp_o.copy()
        current_PS = 0
        new_PSs = []
        sw = tmp.swapped_by_shapeit.iloc[0]
        for i,row in tmp.iterrows():
            if row.swapped_by_shapeit == sw:
                new_PSs.append(current_PS)
            else:
                sw = ~sw
                current_PS = current_PS+1
                new_PSs.append(current_PS)
        return tmp.assign(new_PS = new_PSs)
    # Apply processing
    tmp2 = shpm[shpm.PS.isin(strange_PS)].groupby('PS')['Chromosome POS REF ALT PS_copy swapped_by_shapeit'.split()].apply(process_df).reset_index()
    output = tmp2.sort_values(by='level_1')
    output = output.loc[:,'Chromosome POS REF ALT new_PS'.split()]
    whp = pd.read_csv(whatshap_vcf, sep='\t', comment='#', header=None)

    # Rename columns
    columns = 'Chromosome POS ID REF ALT Score FILTER INFO FORMAT Sample'.split()
    whp = whp.rename(columns={i: x for i, x in enumerate(columns)})
    whp = whp.assign(Chromosome = whp.Chromosome.astype(str))
    whp = whp.assign(PS = [x.split(':')[-1] for x in whp.Sample])
    whp = whp.assign(PS_copy = whp.PS)
    whp = whp.merge(output, on = 'Chromosome POS REF ALT'.split(), how='left')

    def process_tst(starts, pos, vals):
        new_vals = []
        if starts[0] !=0:
            starts = [0]+ starts
            vals = [0.0]+vals
        new_val = pos[0]
        for i in range(len(pos)):
            if i==len(starts):
                starts = starts+[i]
                new_vals.append(new_val)
            else:
                if starts[i] == i:
                    if (i>0):
                        if vals[i]>vals[i-1]:
                            new_val = pos[i]
                    new_vals.append(new_val)
                else:
                    #find the val
                    position = pos[i]
                    left = pos[i-1]
                    right = pos[starts[i]]
                    if vals[i-1] == vals[i]:
                        chosen_val = vals[i]
                    else:
                        if position - left <= right-position:
                            chosen_val = vals[i-1]
                        else:
                            chosen_val = vals[i]
                            new_val = pos[i]
                    vals =  vals[:i]+[chosen_val]+vals[i:]
                    starts = starts[:i]+[i]+starts[i:]
                    new_vals.append(new_val)
        return new_vals

    def process_df_2(tst):
        starts = [i for i in range(len(tst)) if not pd.isna(tst.new_PS.iloc[i])]
        if len(starts)==0:
            return tst.assign(newer_PS = tst.PS_copy)
        pos = list(tst.POS)
        vals = [i for i in tst.new_PS if not pd.isna(i)]
        return tst.assign(newer_PS = process_tst(starts, pos, vals))
    whp = whp.groupby('PS').apply(process_df_2, include_groups=False)
    whp= whp.drop(columns = 'PS_copy new_PS'.split()).reset_index().drop(columns='PS')
    whp = whp.assign(Sample = [':'.join(x.Sample.split(':')[:-1]+[str(int(x.newer_PS))]) for i,x in whp.iterrows()])
    whp = whp.drop(columns = 'newer_PS')
    whp.sort_values(by='level_1').drop(columns='level_1').to_csv(output_file, sep='\t', index=False, quoting=csv.QUOTE_NONE)

    print(f"Output saved to {output_file}")
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF files")
    parser.add_argument("shapeit_vcf", help="Path to the Shapeit VCF file")
    parser.add_argument("whatshap_vcf", help="Path to the Whatshap VCF file")
    parser.add_argument("output_file", help="Path to save the output TSV file")

    args = parser.parse_args()

    process_vcf(args.shapeit_vcf, args.whatshap_vcf, args.output_file)
