from Bio import SeqIO
from Bio.Seq import Seq

dnatypelist = ['1. The Standard Code',
            '2. The Vertebrate Mitochondrial Code',
            '3. The Yeast Mitochondrial Code',
            '4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code',
            '5. The Invertebrate Mitochondrial Code',
            '6. The Ciliate, Dasycladacean and Hexamita Nuclear Code',
            '9. The Echinoderm and Flatworm Mitochondrial Code',
            '10. The Euplotid Nuclear Code',
            '11. The Bacterial, Archaeal and Plant Plastid Code',
            '12. The Alternative Yeast Nuclear Code',
            '13. The Ascidian Mitochondrial Code',
            '14. The Alternative Flatworm Mitochondrial Code',
            '16. Chlorophycean Mitochondrial Code',
            '21. Trematode Mitochondrial Code',
            '22. Scenedesmus obliquus Mitochondrial Code',
            '23. Thraustochytrium Mitochondrial Code',
            '24. Rhabdopleuridae Mitochondrial Code',
            '25. Candidate Division SR1 and Gracilibacteria Code',
            '26. Pachysolen tannophilus Nuclear Code',
            '27. Karyorelict Nuclear Code',
            '28. Condylostoma Nuclear Code',
            '29. Mesodinium Nuclear Code',
            '30. Peritrich Nuclear Code',
            '31. Blastocrithidia Nuclear Code',
            '33. Cephalodiscidae Mitochondrial UAA-Tyr Code']
startcodons = ['ATG', 'GTG', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'TTA']

def fastaparse(files):
    sequences = []
    fafiles = files.split(";")
    for f in fafiles:
        try:
            for seq_record in SeqIO.parse(f,"fasta"):
                print(seq_record.id)
                pair = (seq_record.id,seq_record.seq)
                sequences.append(pair)
        except:
            print("Something went wrong")

    return sequences

def findorfs(pair, dnatable, min_protein_length, startcod):
    answer = []
    startaalist = []
    for codon in startcod:
        codseq = Seq(codon)
        startaa = codseq.translate(dnatable)
        startaalist.append(startaa)
    startaalist = list(dict.fromkeys(startaalist))
    startaalist = [str(aa) for aa in startaalist]
    seq = pair[1]
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            frame_seq = nuc[frame:]
            trans = str(frame_seq.translate(dnatable))
            aa_end = 0
            if startaalist == []:
                aa_start = 0
            else:
                aa_start = getstartloc(aa_end, trans, startaalist)

            while aa_start < len(trans):
                print(len(trans))
                currentstartcodon = frame_seq[(aa_start * 3):(aa_start * 3 + 3)]
                if currentstartcodon in startcod or startcod == []:
                    pass
                else:
                    aa_end = aa_start+1
                    aa_start = getstartloc(aa_end,trans,startaalist)
                    if aa_start == -1:
                        break
                    else:
                        continue
                aa_end = trans.find("*", aa_start)
                aa_end2 = trans.find("X", aa_start)
                if aa_end2 == -1 and aa_end == -1:
                    break
                if aa_end2 < aa_end:
                    aa_end = aa_end2
                if aa_start == -1:
                    break
                if aa_end == -1:
                    break
                if aa_end - aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame + aa_start * 3 + 1
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3 + 1
                        end = seq_len - frame - aa_start * 3
                    startpos = start - 1
                    answer.append((start, end, strand, seq[startpos:end], trans[aa_start:aa_end]))
                if startaalist == []:
                    aa_start = aa_end+1
                else:
                    aa_start = getstartloc(aa_end, trans, startaalist)
    answer.sort()
    return answer

def getstartloc(aa_end, trans,startaalist):
    startoptions = []
    for aa in startaalist:
        aa_startoption = trans.find(aa, aa_end)
        startoptions.append(aa_startoption)
    aa_start = min(startoptions)
    return aa_start

