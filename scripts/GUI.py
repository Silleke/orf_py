import os
import PySimpleGUI as sg
import functions as fc

if not os.path.exists('../output'):
    os.makedirs('../output')
sg.theme('GreenMono')


# ------ Menu Definition ------ #
menu_def = [['&File', ['&Open', '&Save', 'E&xit', 'Properties']],
            ['&Edit', ['Paste', ['Special', 'Normal', ], 'Undo'], ],
            ['&Help', '&About...'], ]

layout = [  [sg.Menu(menu_def, tearoff=True)],
            [sg.Text('Select one or multiple valid .fasta files')],
            [sg.FilesBrowse(),sg.Input()],
            [sg.Text('-'*170)],
            [sg.Text('Select the correct genetic code for your data \n(These and their numbering relate to those of the NCBI database)')],
            [sg.InputCombo(fc.dnatypelist, "1.The Standard Code",readonly=True, key='dnatype')],
            [sg.Text('-'*170)],
            [sg.Text('Minimum protein length for ORF:'), sg.InputText(100, key='protlen')],
            [sg.Text('-'*170)],
            [sg.Text('Start codons:'), sg.Checkbox('ATG', default=True, key = 'ATG'), sg.Checkbox('GTG', key = 'GTG'), sg.Checkbox('TTG', key = 'TTG'),
            sg.Checkbox('CTG', key = 'CTG'), sg.Checkbox('ATT', key = 'ATT'), sg.Checkbox('ATC', key = 'ATC'), sg.Checkbox('ATA', key = 'ATA'), sg.Checkbox('TTA', key = 'TTA')],
            [sg.Text('-'*170)],
            [sg.Text('Output filename:'), sg.InputText('default', size = (15,1), key = 'filename'), sg.Text('.csv')],
            [sg.Submit(key = 'submit')] ]

window = sg.Window('ORFpy', layout)
# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event in (None, 'Exit'):
        break
    if event == 'submit':
        print(f'you clicked {event}')
        print(f'you chose file {values["Browse"]}')
        fasta = fc.fastaparse(values["Browse"])
        dnatype = values['dnatype']
        min_protein_length = int(values["protlen"])
        dnatable = dnatype.split(".")[0]
        startcod = []
        for codon in fc.startcodons:
            if values[codon] == True:
                startcod.append(codon)
        if fasta:
            file1 = open("../output/%s.csv"%values['filename'], "w+")
            file1.write("Header, ORF length, Start, Stop, Strand, Nucleotide Seq, Protein Seq\n")
            for pair in fasta:
                orf_list = fc.findorfs(pair,dnatable,min_protein_length, startcod)
                for start, end, strand, nucl, pro in orf_list:
                    file1.write(
                        "%s, %i, %i, %i, %i, %s, %s"
                        % (pair[0], len(pro), start, end, strand, nucl, pro)
                    )
                    file1.write("\n")
        else:
            print("Something went wrong")
        file1.close()
        print("file written")
window.Close()
print("end of script")
