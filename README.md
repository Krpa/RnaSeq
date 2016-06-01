# RnaSeq

Mentor: prof.dr.sc. Mile Šikić

Tema: Poravnanje RNA očitanja na poznate gene

Projekt u sklopu završnog rada za preddiplomski studij na [FER](http://www.fer.unizg.hr/)-u.
Skripte su namijenjene kao dodatak već postojećem alatu za poravnanje [GraphMap](https://github.com/isovic/graphmap)
kako bi se omogućilo poravnanje RNA očitanja na poznate gene.

## Upute za korištenje

###transExtract.py
Skripta iz zadane [GTF](http://mblab.wustl.edu/GTF2.html) datoteke generira transkriptom sekvence zadane u 
[FASTA](http://genetics.bwh.harvard.edu/pph/FASTA.html) formatu. Transkriptom se zapisuje kao FASTA datoteka
s više sekvenci.
Skripta očekuje najmanje 3 argumenta: putanju do ulazne GTF datoteke, putanju do ulazne FASTA datoteke i putanju do izlazne
FASTA datoteke. Nakon toga slijede zastavice koje ne moraju biti uključene. 

######Zastavice
* __s__ - uključivanje stranda kod generiranja transkriptoma, skripta će za eksone označene "-" strandom uzimati njihov
reverzni komplement
* __alt__ - uključuje "alternative splicing", skripta će za svaki transkript generirati sve podskupove njegovih eksona (osim
praznog skupa) i za svaki podskup generirati novi transkript. Transkript koji je sastavljen od svih eksona iz originalnog
skupa koji mu pripada neće imati posebnu oznaku, a svi ostali (koji su generirani uzimajući neki pravi podskup) će uz originalno
ime transkripta imati pridodan sufiks "_bitmask A" gdje je bitmask niz duljine broja eksona u originalnom skupu sastavljen od 0 i 1 gdje 1 označava uključivanje eksona,
a 0 isključivanje.

#####Primjer pokretanja
```
python ./RnaSeq/src/transExtract.py ./ulaz/ulaz.gtf ./ulaz/chr_s288c_1.fa ./izlaz/izlaz.fa 
python ./RnaSeq/src/transExtract.py ./ulaz/ulaz.gtf ./ulaz/chr_s288c_1.fa ./izlaz/izlaz_a.fa s alt
```


###transToGenome.py
Skripta iz zadanih [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) i GTF datoteka vraća mapiranje zapisano
u SAM datoteci iz prostora transkriptoma u prostor genoma i zapisuje u datoteku koja je predana kao treći argument. 
Skripta očekuje da je mapiranje rađeno bez uključivanja stranda (da je transExtract.py bio pokrenut bez zastavice __s__).

#####Primjer pokretanja
```
python ./RnaSeq/src/transToGenome.py ./ulaz/sim1_20perc_X5.sam ./ulaz/ulaz.gtf ./izlaz/out.sam
```

## Rad s GraphMap-om
Alat __GraphMap__ se koristi za poravnanje RNA očitanja na transkriptom generiran skriptom __transExtract.py__. 
Dobiveni izlaz se prosljeđuje skripti __transToGenome.py__ koja vraća mapiranje iz prostora transkriptoma u prostor genoma.

#####Primjer pokretanja
```
python ./RnaSeq/src/transExtract.py ./ulaz/ulaz.gtf ./ulaz/chr_s288c_1.fa ./izlaz/izlaz_a.fa s alt
graphmap align -r ./izlaz/izlaz_a.fa -d /ulaz/sim1_20perc_X5.fastq -o ./ulaz/sim1_20perc_X5.sam
python ./RnaSeq/src/transToGenome.py ./ulaz/sim1_20perc_X5.sam ./ulaz/ulaz.gtf ./izlaz/out.sam
```
