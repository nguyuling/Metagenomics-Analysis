## Metagenomics Analysis Pipeline

### Environments and packages
(refer `qc.yml` and `biobakery.yml`)
1. **qc**
    - `fastqc`
    - `multiqc`
    - `cutadapt`
    - `sortmerna`
2. **biobakery**
    - `python=3.10`
    - `metaphlan`
    - `krona`
    - `humann=3.9=py_0`
    - `pandas`
    - `openpyxl`
```bash
cd <ENV_PATH>
    conda create -n <ENV_NAME> -c bioconda -c conda-forge fastqc multiqc cutadapt sortmerna
    conda create -n <ENV_NAME> -c conda-forge -c bioconda -c biobakery -c defaults python=3.10 metaphlan krona humann=3.9 pandas openpyxl
```

### Database
1. **sortmeRNA database **[download-here](https://github.com/sortmerna/sortmerna/tree/master/test_data/refs) — reference rRNA FASTA files curated from SILVA (archaeal, bacterial and eukaryotic) and Rfam (5s and 5.8S rRNA subunits)
2. **MetapPhlAn database** — MetaPhlAn 4 ChocoPhlAn SGB (Species-level Genome Bins) marker database vJun23
```bash
    conda activate <BIOBAKERY_ENV>
    metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307
```
3. **HUMAnN** [The Huttenhower Lab](https://huttenhower.sph.harvard.edu/humann/) — 
```bash — includes chocophlan (nucleotide pangenome database), UniRef90 (database of clustered protein sequences), utility mapping (map UniRef id to other classification systems)
    conda activate <BIOBAKERY_ENV>
    PATH=<ENV_PATH>/lib/python3.10/site-packages/humann/data
    humann_databases --download chocophlan full $PATH --update-config yes
    humann_databases --download uniref uniref90_diamond $PATH --update-config yes
    humann_databases --download utility_mapping full $PATH --update-config yes
```

### Guide to Run the Pipeline
1. Run QC `fastqc` and `multiqc` on all paired-end reads
```bash
    conda activate <QC_ENV>
    bash qc.sh all          # run all unfinished steps, OR
    bash qc.sh allforce     # run all steps regardlessly
```
2. Edit `cutadapt.txt`
    - detected reads pair(s), default flags and their values are written to the `cutadapt.txt`
    - any missing paired-end reads, please add sample_name, r1 and r2 manualy in a new row
    - based on multiqc report, specifies trimming parameter for each reads pair by:
        - editing cutadapt flag on header row (`cutadapt --help` to refer all flags)
        - editing the values of each column (flag) and row (reads pair)
3. Resume QC in terminal [y/n], to perform `cutadapt` and `sortmerna`, pair by pair
```bash
y
```
4. Run bioBakery taxonomic and functional profiling, pair by pair
```bash
    conda activate <BIOBAKERY_ENV>
    bash qc.sh all          # run all unfinished steps, OR
    bash qc.sh allforce     # run all steps regardlessly
```

### Output Files Structure