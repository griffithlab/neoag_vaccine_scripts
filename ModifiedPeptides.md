## Generating the peptide_table.tsv file


The input is a csv file which has the names and base sequences. The name does not have to be unique.
For example, the gene name to indentify each sequence. The sequences should not contian any modifications. 

```bash
Base Sequence name,Base sequence
CUL9,RMLDYYEEISAGDEGEFRQS
CUL9,RVRMLDYYEEISAGDEGEFRQSN
CUL9,RVRMLDYYEEISAGDEGEFR
EXOC4,SVIRTLSTIDDVEDRENEKGR
EXOC4,ISVIRTLSTIDDVEDRENEKGR
EXOC4,LISVIRTLSTIDDVEDRENEKGR
VPS13B,GLRQGLFRLGISLLGAIAGIVD
VPS13B,GEGLRQGLFRLGISLLGAIAG
VPS13B,SLGEGLRQGLFRLGISLLGAI
DYNC1H1,KRFHATISFDTDTGLKQALET
DYNC1H1,GKRFHATISFDTDTGLKQALET
DYNC1H1,KRFHATISFDTDTGLKQAL
MYO9A,FDWIVFRINHALLNSKVLEHNTK
MYO9A,FDWIVFRINHALLNSKVL
MYO9A,SALFDWIVFRINHALLNSKVLEHN
EPG5,KELPLYLWQPSTSEIAVIRDW
```

Using the Docker and executing this command will produce the peptide_table.tsv used by pVAC bind. 
The -n argument is the maximum number of modified peptides and the -m argument is a path to the csv file.

```
docker run -it -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:latest /bin/bash

RunDocker griffithlab/neoang_scripts
python3 /opt/scripts/modify_peptides.py -n 3 -m 038_ModifiedPeptides.csv
```

For example, if you speficify -n 1 then the modified sequences produced will me:

```bash
CUL9.1.n-term-K	KRMLDYYEEISAGDEGEFRQS	K|RMLDYYEEISAGDEGEFRQS
CUL9.1.c-term-K	RMLDYYEEISAGDEGEFRQSK	RMLDYYEEISAGDEGEFRQS|K
CUL9.1.n-term-R	RRMLDYYEEISAGDEGEFRQS	R|RMLDYYEEISAGDEGEFRQS
CUL9.1.c-term-R	RMLDYYEEISAGDEGEFRQSR	RMLDYYEEISAGDEGEFRQS|R
.
.
.
EPG5.n-term-K	KKELPLYLWQPSTSEIAVIRDW	K|KELPLYLWQPSTSEIAVIRDW
EPG5.c-term-K	KELPLYLWQPSTSEIAVIRDWK	KELPLYLWQPSTSEIAVIRDW|K
EPG5.n-term-R	RKELPLYLWQPSTSEIAVIRDW	R|KELPLYLWQPSTSEIAVIRDW
EPG5.c-term-R	KELPLYLWQPSTSEIAVIRDWR	KELPLYLWQPSTSEIAVIRDW|R
```
