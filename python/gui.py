__author__ = 'Ariel'

import logging
import tkinter
# http://biopython.org/DIST/docs/tutorial/Tutorial.html
# http://biopython.org/DIST/docs/api/Bio-module.html
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
Entrez.email = "arielvina@yahoo.es"
# http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
# http://biopython.org/wiki/SeqIO
# http://biopython.org/DIST/docs/api/Bio.SeqRecord.SeqRecord-class.html
# http://biopython.org/DIST/docs/api/Bio.SeqFeature.SeqFeature-class.html
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc16
# from Bio import SeqIO
from Bio import GenBank
# http://tkinter.unpythonic.net/wiki/tkFileDialog
from tkinter import filedialog
from tkinter import scrolledtext
align_file_name = None  # '../alignment/HEV.fas'  # or None to ask first


class Seq_pos:

    def __init__(self, seq_beg=0, seq_end=0):
        self.beg = seq_beg
        self.end = seq_end

    def __str__(self):
        return f"Seq pos: {self.beg}-{self.end}"

    def expand(self, pos):
        self.beg = min(pos.beg, self.beg)
        self.end = max(pos.end, self.end)


class Q_hit_pos:

    def __init__(self, q_pos=Seq_pos(), h_pos=Seq_pos()):
        self.q = q_pos
        self.h = h_pos

    def __str__(self):
        return f"Hit pos: Query {self.q}, Hit{self.h}"

    def adjust_h(self, hit):
        self.h.beg = self.q.beg + (hit.q.beg - hit.h.beg)
        self.h.end = self.h.beg + hit.h.end

    def adjust_q(self, hit):
        self.q.beg = self.h.beg + (hit.h.beg - hit.q.beg)
        self.q.end = self.q.beg + hit.q.end


class App(tkinter.Frame):

    def __init__(self):
        tkinter.Frame.__init__(self, tkinter.Tk(),  width=600, height=600)
        self.master.title('Adding new sequences')
        # http://effbot.org/tkinterbook/pack.htm
        # http://effbot.org/tkinterbook/grid.htm#Tkinter.Grid.grid-method
        # http://infohost.nmt.edu/tcc/help/pubs/tkinter/web/grid.html
        self.grid(sticky=tkinter.NS)

        self.refLen = 0
        self.refSeq = dict()
        self.newSeq = dict()

        self.winfo_toplevel().rowconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        w = 15
        h = 40

        self.ID_original = ID_list(self, "Load original list of ID", w, h)
        self.ID_original.grid(row=0, column=0, sticky=tkinter.NSEW)

        self.ID_add = ID_list(self, "Load ID to add", w, h)
        self.ID_add.grid(row=0, column=1, sticky=tkinter.NSEW)
        tkinter.Button(self, text="BLAST",
                             command=self.blast)                     .grid(row=2, column=1)

        self.ID_unique = ID_list(self, "Load", w, h)
        self.ID_unique.grid(row=0, column=2, sticky=tkinter.NSEW)
        tkinter.Button(self, text="Load BLAST",
                             command=self.load_blast)                .grid(row=2, column=2)
        tkinter.Button(self, text="Seq from GB file",
                             command=self.parseGB)                   .grid(row=2, column=0)
        tkinter.Button(self, text="Filter",
                             command=self.filter)                    .grid(row=3, column=2)
        tkinter.Button(self, text="Parse alignment",
                             command=self.parseAlign)                .grid(row=3, column=0)
        if align_file_name:
            self.parseAlign(align_file_name)

    def parseAlign(self, file_name=None):
        if not file_name:
            file_name = filedialog.askopenfilename( filetypes=(("fasta aligment", "*.fas"), ("All files", "*.*")),
                                                      defaultextension='fas',
                                                      title='Select Master Alignment')
            if not file_name:
                return

        self.refSeq.clear()
        logging.info(file_name)
        with open(file_name) as align_file:
            seq_name=''
            for line in align_file.readlines():
                logging.debug(line)
                if line[0] == '>':
                    seq_name = line[1:].rstrip()
                else:
                    ln = len(line)
                    self.refLen = max(ln, self.refLen)
                    seq_beg = 0
                    seq_end = ln-2
                    while seq_beg<ln :
                        if line[seq_beg] == '-':
                            seq_beg += 1
                        else:
                            break   # todo :  check it is a valid base not line end???
                    while seq_end > seq_beg:
                        if line[seq_end] == '-':
                            seq_end -= 1
                        else:
                            break   # todo :  check it is a valid base not line end???
                    logging.debug('Seq: ' + seq_name + str((seq_beg,seq_end)))
                    self.refSeq[seq_name] = Seq_pos(seq_beg, seq_end)
        self.ID_original.clear()
        self.ID_original.add('\n'.join(self.refSeq.keys()))

    def filter_add(self, add):
        ori  = set([oID.split('.')[0] for oID in self.ID_original.lines()])
        uniq = set([oID.split('.')[0] for oID in self.ID_unique.lines()])
        uniq.update([oID.split('.')[0] for oID in add])
        uniq -= ori
        self.ID_add.clear()
        self.ID_unique.clear()
        for ID in uniq:
            self.ID_unique.add(ID)

    def filter(self):
        self.filter_add(self.ID_add.lines())

    def blast(self, blast_file_name=None):
        """
        Take the set of the ID in the "add" list (center list) and make on-line an NCBI.BLAST.
        Ask for a file to save the results of the BLAST.
        Call self.load_blast_data(blast_file) to parse the results.
        :return:
        """
        IDs = list(set(self.ID_add.lines()))
        logging.debug(' '.join(IDs))
        if not blast_file_name:
            blast_file_name = filedialog.asksaveasfilename(filetypes=(("BLAST", "*.xml"), ("All files", "*.*") ), defaultextension='xml', title='Save the BLAST result in XML format')
        if not blast_file_name:
            return
        # lets limit the number of sequences to BLAST in one pass to NS
        NS = 2
        i = 0
        while i < len(IDs):
            self.master.title('Talking to NCBI. Running BLAST. Be VERY patient ...')
            logging.info('BLAST: ' + ' '.join(IDs[i : i + NS]))
            result_handle = NCBIWWW.qblast("blastn", "nt", '\n'.join(IDs[i:i+NS]))   #, hitlist_size=50, perc_ident=90, threshold=1, alignments=50, filter="HEV", format_type='XML', results_file=blast_file_name )
            self.master.title('Adding new sequences...')
            print('... Returned from NCBI BLAST')
            # http://tkinter.unpythonic.net/wiki/tkFileDialog
            with open(blast_file_name+'-'+str(i), mode='w') as blast_file:
                blast_file.write(result_handle.read())
            result_handle.close()
            # self.load_blast_data(result_handle)
            with open(blast_file_name+'-'+str(i), mode='r') as blast_file:
                self.load_blast_data(blast_file)
            i += NS

    def blast_seq(self, seq, entrez_query="(none)", blast_file_name=None):
        if not blast_file_name:
            blast_file_name = filedialog.asksaveasfilename(filetypes=(("BLAST", "*.xml"), ("All files", "*.*") ),
                                                           defaultextension='xml', title='Save the BLAST result in XML format')
        if not blast_file_name:
            return
        self.master.title('Talking to NCBI. Running BLAST. Be VERY patient ...')
        logging.info('BLAST: ' + seq)
        result_handle = NCBIWWW.qblast(program="blastn",
                                       database="nt",
                                       sequence=seq,
                                       entrez_query=entrez_query,
                                       descriptions=5000,
                                       alignments=5000,
                                       hitlist_size=5000
                                       )  # , hitlist_size=50, perc_ident=90, threshold=1, alignments=50, filter="HEV",
                                          # format_type='XML', results_file=blast_file_name )

        self.master.title('Adding new sequences...')
        print('... Returned from NCBI BLAST')
        with open(blast_file_name, mode='w') as blast_file:
            blast_file.write(result_handle.read())
        result_handle.close()
        # self.load_blast_data(result_handle)
        with open(blast_file_name, mode='r') as blast_file:
            self.load_blast_data(blast_file)


    def load_blast(self):
        with filedialog.askopenfile(filetypes=(("BLAST (xml)", "*.xml"), ("All files", "*.*") ),
                                    title='Load a BLAST result in XML format') as blast_file:
            logging.info(f'Proccesing BLASt file: {blast_file.name}')
            self.load_blast_data(blast_file)


    def load_blast_data(self, blast_data):
        """
        Uses NCBIXML.parse to parse the result of the BLAST and identify each sequence there:
         Look for each sequence in the self.refSeq dictionary to know if it is a reference sequence (the ones added with
         parse alignment) and add it to self.newSeq if not a ref.
        :param blast_data: an open XML BLAST file.
        :return:
        """
        IDs = set()
        blast_records = NCBIXML.parse(blast_data)
        logging.info('proccesing BLASt file: parsed')
        q_is_ref = False
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc91
        for blast_record in blast_records:
            try:
                qID = blast_record.query_id.split('|')[3].split('.')[0]
            except:
                qID = blast_record.query
            q_is_ref = qID in self.refSeq.keys()
            align_pos = Q_hit_pos()
            if q_is_ref:
                align_pos.q = self.refSeq[qID]
            h_is_ref = False
            logging.debug('Query: ' + blast_record.query_id + ' : ' + qID + '. Is a Reference: ' + str(q_is_ref))
            for alignment in blast_record.alignments:
                h_is_ref = alignment.accession in self.refSeq.keys()
                logging.debug('Hit: ' + str(alignment.accession) + '. Is a Reference: ' + str(h_is_ref))
                IDs.add(alignment.accession)           # alignment.title.split('|')[3].split('.')[0])
                if q_is_ref and h_is_ref:
                    continue
                h = self.hit_positions(alignment)
                logging.debug(h)
                logging.debug(align_pos)
                if q_is_ref:
                    if not h_is_ref:
                        align_pos.adjust_h(h)
                        if alignment.accession in self.newSeq:
                            self.newSeq[alignment.accession].expand(align_pos.h)
                        else:
                            self.newSeq[alignment.accession] = align_pos.h
                else:
                    if h_is_ref:
                        align_pos.h = self.refSeq[alignment.accession]
                        align_pos.adjust_q(h)
                        if qID in self.newSeq:
                            self.newSeq[qID].expand(align_pos.q)
                        else:
                            self.newSeq[qID] = align_pos.q


                logging.debug(align_pos)

        self.filter_add(IDs)

    def load_primer_blast(self):
        with filedialog.askopenfile(filetypes=(("BLAST Results", "*.html"), ("All files", "*.*") ),
                                    title='Load a Primer-BLAST result in HTML format') as primer_blast_file:
            logging.info(f'Proccesing Primer-BLAST file: {primer_blast_file.name}')
            self.load_primer_blast_data(primer_blast_file)

    def load_primer_blast_data(self, primer_blast_file):
        IDs = set()
        seq_name = ''
        new_entrez = False
        while not new_entrez:  # skip html headers
            line = primer_blast_file.readline()
            logging.debug(f'Parsing line: {line}')
            if not line: return
            new_entrez = line.split('target="new_entrez">')
            wwwt = new_entrez.pop(0)

        while True:
            while not new_entrez:  # skip html headers
                line = primer_blast_file.readline()                               # new_entrez ? or just <pre>?
                if not line: return
                if not "<pre>" == line[:-1]:
                    new_entrez = line.split('target="new_entrez">')
                    wwwt = new_entrez.pop(0)
                    continue      # <pre>
                logging.debug(f'Parsing line: {line}')
                product_lenght = int(primer_blast_file.readline().split('=')[1])   # product length = 102
                line = primer_blast_file.readline()                                # Forward primer  1       GCCTTCCAGACCATGCTC  18
                _, left, pattern, right = primer_blast_file.readline().split()     # Template        202309  ..................  202326
                line = primer_blast_file.readline()                                #
                line = primer_blast_file.readline()                                # Reverse primer  1       AGTGCGGAGGTCATTTGC  18
                _, rleft, rpattern, rright = primer_blast_file.readline().split()  # Template        202410  ..................  202393
                line = primer_blast_file.readline()                                #
                line = primer_blast_file.readline()                                # </pre>

                print (f"{acc}\tAmplify\t{product_lenght}\tnt\tfw:\t{pattern}\t{left}\t{right}\trv:\t{rpattern}\t{rleft}\t{rright}\t{desc[:-1]}")

            logging.debug(f'Parsing new_entrez: {new_entrez}')
            acc, desc = new_entrez[0].split('</a> ')
            www = wwwt
            new_entrez = False

        # print('proccesing BLASt file: parsed')

    def hit_positions(self, alignments):
        #assert (isinstance(alignments, NCBIXML.alignment))
        pos=[Q_hit_pos(Seq_pos(hit.query_start, hit.query_end), Seq_pos(hit.sbjct_start, hit.sbjct_end)) for hit in alignments.hsps]
        h = pos[0]
        if len(pos) == 1:
            return h
        for p in pos[1:]:
            if  h.q.beg < p.q.beg:
                h.q.beg = p.q.beg
                h.h.beg = p.h.beg
            if  h.q.end > p.q.end :
                h.q.end = p.q.end
                h.h.end = p.h.end
        return h

    def get_seq_GB(self, IDs):
        '''
        Get the sequences corresponding to the list of IDs from GenBank online,
        save the sequences in a flat file and call the parser
        Horror !! Surprise ! Biopython don't support sequence GenBank xml files , only flat.
        :return:
        '''

        logging.debug (IDs)
        seq_flat_file_name = filedialog.asksaveasfilename(filetypes=(("Seq flat GB", "*.gb"), ("All files", "*.*")),
                                                          defaultextension='gb',
                                                          title='Save the GenBank sequences in flat format')

        NS = 100
        i = 0
        IDs = IDs.split(',')
        with open(seq_flat_file_name, mode='w') as seq_file:
            while i < len(IDs):
                #ID =
                logging.debug('GenBank: ' + str(IDs[i: i + NS]))
                self.master.title('Talking to the NCBI. Getting sequences. Be VERY patient ...')
                seq_handle = Entrez.efetch(db="nuccore",
                                           id=IDs[i: i + NS],
                                           rettype="gb",
                                           retmode="text" )
                # Entrez.efetch(db="nucleotide", id="57240072", rettype="gb", retmode="text")
                i += NS
                self.master.title('Adding new sequences...')
                print('returned')

                for line in seq_handle:
                    seq_file.write(line)
                seq_handle.close()

        self.parseGBfile(seq_flat_file_name)

    def parseGB(self):
        '''
        Load and parse a GenBank sequence flat file
        :return:
        '''
        seq_flat_file_name = filedialog.askopenfilename(filetypes=(("Seq flat GB", "*.gb"), ("All files", "*.*") ),
                                                        title='Parse the GenBank sequences in flat format')
        self.parseGBfile(seq_flat_file_name)

    def parseGBfile(self, seq_flat_file_name):
        # https://docs.python.org/3.4/library/pathlib.html#module-pathlib
        # https://docs.python.org/3.4/library/os.path.html#module-os.path
        fasta_file_name = seq_flat_file_name.replace('.gb', '')+'.fasta'
        csv_file_name = seq_flat_file_name.replace('.gb', '')+'.csv'
        sep = ';'
        el = '\n'
        seq = []
        with open(csv_file_name, 'w') as csv:
            with open(seq_flat_file_name) as seq_flat_file:
                # http://biopython.org/DIST/docs/api/Bio.GenBank-module.html#parse
                for record in GenBank.parse(seq_flat_file):     #, "genbank"
                    #  http://biopython.org/DIST/docs/api/Bio.GenBank.Record-module.html
                    #  http://biopython.org/DIST/docs/api/Bio.GenBank.Record.Record-class.html
                    #  http://biopython.org/DIST/docs/api/Bio.GenBank.Scanner-pysrc.html#GenBankScanner._feed_header_lines
                    if record.locus in self.refSeq:
                        continue
                    beg = 0
                    end = 0
                    if record.locus in self.newSeq:
                        beg, end = self.newSeq[record.locus]
                        record.sequence = '-' * beg + record.sequence
                        record.sequence = record.sequence + '-' * (self.refLen-len(record.sequence) )

                    sq = '>' + record.locus + el + record.sequence + el
                    # fasta.write('>' + record.locus + el + record.sequence +el)   # record.accession[0]  ??
                    seq.append((sq, beg))
                    csv.write(record.locus + sep)               # MEGA name:(A)
                    csv.write('no'         + sep)               # Tab-Pub:  (B)

                    strain  = ''
                    isolate = ''
                    host    = ''
                    country = ''
                    region  = ''
                    collection_date = ''
                    source  = ''
                    genotype = ''
                    #  http://biopython.org/DIST/docs/api/Bio.GenBank.Record-pysrc.html#Feature
                    for feature in record.features:
                        if feature.key == 'source':
                            # http://biopython.org/DIST/docs/api/Bio.GenBank.Record.Qualifier-class.html
                            for q in feature.qualifiers:
                                if q.key == '/strain=':
                                    strain = q.value[1:-1].strip()
                                elif q.key == '/isolate=':
                                    isolate = q.value[1:-1].strip()
                                elif q.key == '/country=':
                                    country = q.value[1:-1].split(':')
                                    if len(country) > 1:
                                        region = country[1].strip() # ok?
                                    country = country[0].strip()
                                elif q.key == '/collection_date=':
                                    collection_date = q.value[1:-1].strip()
                                elif q.key == '/source=' or q.key == '/isolation_source=':
                                    source = q.value[1:-1].strip()
                                elif q.key == '/note=':
                                    val = q.value[1:-1]
                                    if val.startswith('genotype: '):
                                        genotype = val[9:].strip()
                    des = record.definition
                    if 'isolate' in des:
                        iso = des.split('isolate')[1]
                        if iso[0] == ':':
                            iso = iso[1:].strip()
                        iso = iso.split()[0].strip()
                        if iso[-1] == '.':
                            iso = iso[:-1].strip()
                        if isolate == '':
                            isolate = iso
                        else:
                            if isolate != iso:
                                isolate = isolate + ' or ' + iso
                    if 'strain' in des:
                        st = des.split('strain')[1]
                        if st[0] == ':':
                            st = st[1:].strip()
                        st = st.split()[0].strip()
                        if st[-1] == '.':
                            st = st[:-1].strip()
                        if strain == '':
                            strain = st
                        else:
                            if strain != st:
                                strain = strain + ' or ' + st


                    csv.write(genotype + sep)            #           ( C )
                    csv.write(sep+sep)                   #           ( C D E)
                    csv.write(strain  +sep)              # Strain name: (F)
                    csv.write(isolate +sep)              # isolate:   (G )
                    csv.write(country +sep +sep +sep)    # country:   (H I J)
                    csv.write(sep + region + sep)        # country:   (KL)
                    csv.write(host +sep)                 # host:      (M)
                    csv.write(source +sep)               # source:    (N)
                    csv.write(collection_date +sep +sep) # year !!! parse !! (OPQ)
                    csv.write(sep + sep)                 #            (RS)
                    csv.write(str(len(record.sequence)) + sep)# Length     (RS)
                    csv.write(record.definition + sep)

                    csv.write(el)

        with open(fasta_file_name, 'w') as fasta:
            seq.sort(key=lambda s : s[1])
            for s in seq:
                fasta.write(s[0])


class ID_list(tkinter.Frame):
    def __init__(self, root, load_titel, width=15, height=40):
        tkinter.Frame.__init__(self, root)
        self.grid(sticky=tkinter.NS)
        self.rowconfigure(1, weight=1)

        tkinter.Button(self, text=load_titel,
                             command=self.load)                      .grid(row=0, column=0, columnspan=3)
        self.txt_list = scrolledtext.ScrolledText(self,  width=width, height=height)
        self.txt_list                                                .grid(row=1, column=0, sticky=tkinter.NSEW, padx=2, pady=2, columnspan=3)
        tkinter.Button(self, text="clear",
                             command=self.clear)                     .grid(row=2, column=0)
        tkinter.Button(self, text="Save",
                             command=self.save)                      .grid(row=2, column=1)
        tkinter.Button(self, text="Get",
                             command=self.get)                       .grid(row=2, column=2 )

    def add(self, ID):
        self.txt_list.insert(tkinter.END, ID+'\n')

    def clear(self):
        self.txt_list.delete(1.0, tkinter.END)

    def load(self):
        with filedialog.askopenfile(filetypes=(("TXT", "*.txt"), ("All files", "*.*") ),
                                        title='Load a ID list in txt format'             ) as ID_file:
            self.add(ID_file.read())
            #for ID in ID_file:
            #    self.add(ID)

    def save(self):
        with filedialog.asksaveasfile(mode='w',
                                      filetypes=(("TXT", "*.txt"), ("All files", "*.*") ),
                                      defaultextension='txt',
                                      title='Save the ID list in txt format') as ID_file:
            ID_file.write('\n'.join(self.lines()))

    def lines(self):
        return self.txt_list.get('1.0',tkinter.END).splitlines()

    def get(self):
        self.master.get_seq_GB(', '.join(self.lines()))


if __name__=='__main__':
    App().mainloop()
