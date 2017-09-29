print('tk...')
from tkinter import filedialog
print('sqlite3...')
import sqlite3
print('Bio...')
from Bio import SeqIO
print('openpyxl...')
import openpyxl
# import tkinter


#class RefSchema


class CreateTaxa:
    def __init__(self, db, kingdom, root_taxa, NCBI_TaxID  = None):
        self.db=db
        self.c = db.cursor()
        self.kingdom = kingdom
        self._root_rank('superkingdom', kingdom)
        self._root_taxa(root_taxa, root_taxa, NCBI_TaxID)

    def _root_rank (self, name, kingdom):        #  todo: NCBI ? add kingdom ! determine type of rank.
        self.c.execute("INSERT INTO taxa_rank (Name, kingdom )"
                       "               VALUES (?   , ?       )",
                       (                       name, kingdom )   )
        self.r_rank = self.c.lastrowid
        return self.r_rank

    def _root_taxa (self, Name ,  vulgar, NCBI_TaxID):
        self.c.execute("INSERT INTO taxa      (Name , vulgar, Id_rank    , NCBI_TaxID  )"
                       "               VALUES ( ?   , ?     , ?          , ?           )",
                                              (Name , vulgar, self.r_rank, NCBI_TaxID  ) )
        self.r_taxa = self.c.lastrowid
        self.c.execute("INSERT INTO taxa_parents (Id_taxa     , parent     , Id_rank      )"
                       "                  VALUES ( ?          , ?          , ?            )",
                                                 (self.r_taxa , self.r_taxa, self.r_rank  ) )
        return self.r_taxa

    def rank (self, name, parent_rank):
        self.c.execute("INSERT INTO taxa_rank (Name, kingdom     , parent      )"
                       "               VALUES (?   , ?           , ?           )",
                       (                       name, self.kingdom, parent_rank )   )
        return self.c.lastrowid

    def taxa (self, Name, vulgar, rank    , parent_taxa,  NCBI_TaxID  = None):
        self.c.execute("INSERT INTO taxa (Name,  vulgar, Id_rank, parent     , NCBI_TaxID )"
                       "          VALUES (?   , ?      , ?      , ?          , ?          )",
                       (                  Name, vulgar , rank   , parent_taxa, NCBI_TaxID )   )
        Id_taxa=self.c.lastrowid
        self.c.execute("INSERT INTO taxa_parents (Id_taxa, parent, Id_rank      )"
                       "                  VALUES ( ?     , ?     , ?            )",
                                                 (Id_taxa , Id_taxa, rank       ) )
        self.c.execute("INSERT INTO taxa_parents (Id_taxa, parent, Id_rank      )"
                       "                  SELECT  ?      , parent, Id_rank       "
                       "FROM taxa_parents WHERE Id_taxa=?",     (Id_taxa , parent_taxa  ) )
        return Id_taxa

def create(newly) -> sqlite3.Connection:
    db = sqlite3.connect("../data/temp/seq.db")
    if newly:
       read_create(db)
       print('Adding default taxas...')
       add_def_taxa(db)
       print('Adding reference schemes...')
       add_ref_schema(db)
       db.commit()

    return db

def read_create(db):
    with open("create_seq.sql") as dbcreate:
        sql_create = dbcreate.read()
    c = db.cursor()
    c.executescript(sql_create)

def read_country_codes(db):
    with open("country_codes.sql") as dbcreate:
        sql_create = dbcreate.read()
    c = db.cursor()
    c.executescript(sql_create)
    db.commit()


def add_ref_schema(db):
    c = db.cursor()
    for sch in [('Lu'  , 'Lu, Li, 2006'         ),
                ('VR'  , 'Vina-Rodriguez, 2015' ),
                ('ICVT', 'ICVT, 2016'           )   ] :
      c.execute("INSERT INTO ref_schema (schema, name) VALUES (?,     ?)", sch)


def add_def_taxa(db):
    ct= CreateTaxa(db, 'viruses', 'Viridae', '10239')

    r   = ct.rank('no rank', ct.r_rank)
    t   = ct.taxa('ssRNA viruses', 'viruses', r, ct.r_taxa, '439488')

    r   = ct.rank('', r)
    t   = ct.taxa('ssRNA positive-strand viruses, no DNA stage', 'viruses', r, t, '35278')

    r   = ct.rank('family', r)
    t   = ct.taxa('Hepeviridae', 'HEV'  , r, t, '291484')

    rGenus  = ct.rank('genus', r)
    tOrth  = ct.taxa('Orthohepevirus', 'Orthohepevirus'  , rGenus, t, '1678141')
    tPisci = ct.taxa('Piscihepevirus', 'Piscihepevirus'  , rGenus, t, '1678142')

    rSpecie  = ct.rank('species', rGenus)
    tOrthSpcA= ct.taxa('Orthohepevirus A', 'Orthohepevirus A'  , rSpecie, tOrth, '1678143')
    tOrthSpcB= ct.taxa('Orthohepevirus B', 'Orthohepevirus B'  , rSpecie, tOrth, '1678144') # NCBI_ID tentative !!
    tOrthSpcC= ct.taxa('Orthohepevirus C', 'Orthohepevirus C'  , rSpecie, tOrth, '1678145')
    tOrthSpcD= ct.taxa('Orthohepevirus D', 'Orthohepevirus D'  , rSpecie, tOrth, '1678146')

    tPisciSpcA=ct.taxa('Piscihepevirus A', 'Piscihepevirus A'  , rSpecie, tPisci,'1678146')

    rGenotype = ct.rank('genotype', rSpecie)
    g1  = ct.taxa('1', 'HEV-g1'  , rGenotype, tOrthSpcA, '185579')
    g2  = ct.taxa('2', 'HEV-g2'  , rGenotype, tOrthSpcA )
    g3  = ct.taxa('3', 'HEV-g3'  , rGenotype, tOrthSpcA, '509628')
    g4  = ct.taxa('4', 'HEV-g4'  , rGenotype, tOrthSpcA, '185580')
    g5  = ct.taxa('5', 'HEV-g5'  , rGenotype, tOrthSpcA )
    g6  = ct.taxa('6', 'HEV-g6'  , rGenotype, tOrthSpcA )
    g7  = ct.taxa('7', 'HEV-g7'  , rGenotype, tOrthSpcA )

    gC1 = ct.taxa('C1', 'HEV-C1' , rGenotype, tOrthSpcC )
    gC2 = ct.taxa('C2', 'HEV-C1' , rGenotype, tOrthSpcC )

    rmc  = ct.rank('major clade', rGenotype)
    maI  = ct.taxa('I'  , 'HEV-g3-I'     , rmc, g3)
    maII = ct.taxa('II' , 'HEV-g3-II'    , rmc, g3)
    Rab  = ct.taxa('3ra', 'HEV-g3-rabbit', rmc, g3)

    rgr  = ct.rank('group', rmc)
    grchi  = ct.taxa('3chi'  , 'HEV-g3chi'     , rgr, maI)
    grjab  = ct.taxa('3jab'  , 'HEV-g3jab'     , rgr, maI)
    grfeg  = ct.taxa('3feg'  , 'HEV-g3feg'     , rgr, maII)

    rsubt = ct.rank('subtype', rgr)
    
    g1a   = ct.taxa('1a', 'HEV-g1a'  , rsubt, g1)
    g1b   = ct.taxa('1b', 'HEV-g1b'  , rsubt, g1)
    g1c   = ct.taxa('1c', 'HEV-g1c'  , rsubt, g1)
    g1d   = ct.taxa('1d', 'HEV-g1d'  , rsubt, g1)
    g1e   = ct.taxa('1e', 'HEV-g1e'  , rsubt, g1)
    g1f   = ct.taxa('1f', 'HEV-g1f'  , rsubt, g1)
    g1g   = ct.taxa('1g', 'HEV-g1g'  , rsubt, g1)
    g1h   = ct.taxa('1h', 'HEV-g1h'  , rsubt, g1)
    g1i   = ct.taxa('1i', 'HEV-g1i'  , rsubt, g1)
    g1j   = ct.taxa('1j', 'HEV-g1j'  , rsubt, g1)
    g1k   = ct.taxa('1k', 'HEV-g1k'  , rsubt, g1)

    g3a   = ct.taxa('3a', 'HEV-g3a'  , rsubt, grjab)
    g3b   = ct.taxa('3b', 'HEV-g3b'  , rsubt, grjab)
    g3c   = ct.taxa('3c', 'HEV-g3c'  , rsubt, grchi)
    g3d   = ct.taxa('3d', 'HEV-g3d'  , rsubt, g3)
    g3e   = ct.taxa('3e', 'HEV-g3e'  , rsubt, grfeg)
    g3ef  = ct.taxa('3ef', 'HEV-g3ef', rsubt, grfeg)
    g3f   = ct.taxa('3f', 'HEV-g3f'  , rsubt, grfeg)
    g3g   = ct.taxa('3g', 'HEV-g3g'  , rsubt, grfeg)
    g3h   = ct.taxa('3h', 'HEV-g3h'  , rsubt, grchi)
    g3i   = ct.taxa('3i', 'HEV-g3i'  , rsubt, grchi)
    g3j   = ct.taxa('3j', 'HEV-g3j'  , rsubt, grjab)
    g3k   = ct.taxa('3k', 'HEV-g3k'  , rsubt, g3)
    g3l   = ct.taxa('3l', 'HEV-g3l'  , rsubt, g3)

    g4a   = ct.taxa('4a', 'HEV-g4a'  , rsubt, g4)
    g4b   = ct.taxa('4b', 'HEV-g4b'  , rsubt, g4)
    g4c   = ct.taxa('4c', 'HEV-g4c'  , rsubt, g4)
    g4d   = ct.taxa('4d', 'HEV-g4d'  , rsubt, g4)
    g4e   = ct.taxa('4e', 'HEV-g4e'  , rsubt, g4)
    g4f   = ct.taxa('4f', 'HEV-g4f'  , rsubt, g4)
    g4g   = ct.taxa('4g', 'HEV-g4g'  , rsubt, g4)
    g4h   = ct.taxa('4h', 'HEV-g4h'  , rsubt, g4)
    g4i   = ct.taxa('4i', 'HEV-g4i'  , rsubt, g4)
    g4j   = ct.taxa('4j', 'HEV-g4j'  , rsubt, g4)
    g4k   = ct.taxa('4k', 'HEV-g4k'  , rsubt, g4)

    db.commit()

def build_ref_pos(Seq, beg, end, Al_len):
    sr=0
    ref = [sr]*beg
    for b in Seq:
        if (b != '-'): sr+=1
        ref.append(sr)
    ref += [sr]*(Al_len - end-1)
    return ref

def parse_full_fasta_Align(db, ref_seq = None, file_name=None):
    """Will parse an alignment and insert it into the tables:
        files: the file path ,
        align:
        seq: because assume this sequences are all complete original sequences (but may be a partial sequence from
             some isolate ! )

        """
    if not file_name:
        file_name = filedialog.askopenfilename(filetypes=(("fasta aligment", "*.fas"), ("All files", "*.*")),
                                               defaultextension='fas',
                                               title='Select Master Alignment')
        if not file_name:
            return

    # self.refSeq.clear()
    print(file_name)

    c = db.cursor()
    c.execute("INSERT INTO seq_file (path, format) VALUES (?, 'fasta')", (file_name,))

    Id_file = c.lastrowid
    c.execute("INSERT INTO align (Id_file, Name,      Ref     ) "
              "           VALUES (?,       ?,         ?       )",
                                 (Id_file, file_name, ref_seq )      )

    Id_align = c.lastrowid
    max_len = 0
    for seq_record in SeqIO.parse(file_name, "fasta"):
        # print(seq_record.id, len(seq_record) )
        seq = str(seq_record.seq)
        if ref_seq is None:      # set first seq as reference
            ref_seq = seq_record.id
            al_ref_seq = seq
        else:
            if ref_seq == seq_record.id:
                al_ref_seq = seq
        ln = len( seq)
        if max_len < ln: max_len = ln
        seq_beg = 0
        seq_end = ln - 1
        while seq_beg < ln:
            if  seq[seq_beg] == '-':
                seq_beg += 1
            else:
                break
        while seq_end > seq_beg:
            if seq[seq_end] == '-':
                seq_end -= 1
            else:
                break

        seq = str( seq[seq_beg: seq_end + 1])
        exp_seq = seq.replace('-','')  #''.join([base for base in seq if base != '-'])

        c.execute("INSERT INTO seq (Name,               Seq,     Len  ) "
                  "         VALUES (?,                  ?,       ?    )",
                                   (str(seq_record.id), exp_seq, len(exp_seq))    )
        Id_part = c.lastrowid

        c.execute("INSERT INTO aligned_seq (Id_align, Id_part, Seq,      pbeg,     pend  ) "
                  "                 VALUES (?,        ?,       ?,        ?,       ?    )",
                                           (Id_align, Id_part, seq, seq_beg, seq_end )    )

    c.execute("UPDATE align SET Al_len = ?, Ref=?    WHERE Id_align=?",
                            (   max_len   , ref_seq,        Id_align ) )
    db.commit()
    return Id_align , build_ref_pos(al_ref_seq,0,len(al_ref_seq)-1, max_len)

def ref_pos(sdb, ID_align, seq_name=None):
    c = sdb.cursor()
    c.execute("SELECT Al_len, Ref FROM align WHERE Id_align=? ", (ID_align,  ))
    Al_len, Ref = c.fetchone()
    if not seq_name: seq_name = Ref
    c.execute("SELECT aligned_seq.Seq, pbeg, pend FROM aligned_seq, Seq "
              "ON Id_part=Id_seq WHERE Name=? AND Id_align=?",
                                     (seq_name,     ID_align ))
    Seq, beg, end = c.fetchone()
    return build_ref_pos(Seq, beg, end, Al_len)

def abnormal_row(c, row):
    MEGA_name = row[0].value          # 'A' - MEGA name
    subtype   = row[3].value          # 'D' - subtype
    Str_name  = row[5].value          # 'F ' 5 - Str.name
    if not subtype: subtype   = row[2].value          # 'C' - genotype

    c.execute("SELECT Id_taxa FROM taxa WHERE taxa.Name=?", (subtype, ))
    Id_taxa = c.fetchone()
    Id_taxa = Id_taxa[0] if Id_taxa else Id_taxa

    c.execute("SELECT Id_seq FROM seq WHERE seq.Name=? ", ( MEGA_name,))
    Id_seq = c.fetchone()
    Id_seq = Id_seq[0] if Id_seq else Id_seq

    c.execute("SELECT Id_algseq FROM aligned_seq WHERE Id_part=?", (Id_seq,))
    Id_algseq= c.fetchone()
    Id_algseq = Id_algseq[0] if Id_algseq else Id_algseq

    if Id_algseq is not None:
        c.execute("INSERT INTO classified_seq (Id_taxa, Id_algseq) VALUES (?,?) "
                                            , (Id_taxa, Id_algseq)                 )
    else:
        c.execute("INSERT INTO pending_seq (Id_taxa, Name,      Id_seq) VALUES (?,?,?)",
                                           (Id_taxa, MEGA_name, Id_seq))

    print("Abnormal row !!!!! ", MEGA_name, subtype, Str_name,"-------> Taxa:{0}, Alseq:{1}, Seq:{2}".format( Id_taxa, Id_algseq, Id_seq))
    if Id_algseq : return False
    return True

def parse_row(db,row, c):
    success = True
    MEGA_name = row[c['MEGA name'   ]].value
    genotype  = row[c['genotype'    ]].value
    subtype   = row[c['subtype'     ]].value
    group     = row[c['group'       ]].value
    Str_name  = row[c['Str.name'    ]].value
    Isolate   = row[c['Isolate'     ]].value
    #Country  = row[c['Country'     ]].value
    Country_cod3=row[c['Country cod']].value
    Region     =row[c['region'      ]].value         # 'H '        todo: deduced
    Region_full=row[c['region full' ]].value
    Host       =row[c['Host'        ]].value
    Source     =row[c['Source'      ]].value
    Year       =row[c['Y'           ]].value
    Month      =row[c['M'           ]].value
    Day        =row[c['D'           ]].value
    Institut   =row[c['Inst'        ]].value
    Reference  =row[c['reference'   ]].value
    Lu_Li      =row[c['Lu, Li'      ]].value
    CG         =row[c['CG'          ]].value

    if not subtype: subtype   = group
    if not subtype: subtype   = genotype

    if not Isolate: Isolate   = Str_name



    c = sdb.cursor()
    c.execute("SELECT Id_strain FROM strain WHERE Name=?", (Str_name, ))
    Id_strain = c.fetchone()
    if Id_strain:
        # print('Existing Strain:', Str_name, Id_strain)
        Id_strain = Id_strain[0]
    else:
        # print('New Strain:', Str_name, Id_strain)
        c.execute("INSERT INTO strain (Name) VALUES (?) ", (Str_name,))
        Id_strain = c.lastrowid
        # print('New Strain ID:', Id_strain)

    c.execute("SELECT Id_taxa FROM taxa WHERE taxa.Name=?", (subtype, ))   # ?? Name UNIQUE ??
    Id_taxa = c.fetchone()
    Id_taxa = Id_taxa[0] if Id_taxa else Id_taxa

    c.execute("SELECT Id_seq FROM seq WHERE seq.Name=? ", ( MEGA_name,))   # ?? Name UNIQUE ??
    Id_seq = c.fetchone()
    Id_seq = Id_seq[0] if Id_seq else Id_seq

    # revise this. Is general??
    c.execute("SELECT Id_algseq FROM aligned_seq WHERE Id_part=?", (Id_seq,))
    Id_algseq= c.fetchone()
    Id_algseq = Id_algseq[0] if Id_algseq else Id_algseq

    # todo: parse location. Is unique?
    c.execute("INSERT INTO isolate (Name   , Id_strain, year, month, day, host, source, institution, country_iso3, region, region_full ) "
              "             VALUES (?      , ?        , ?   , ?    , ?  , ?   , ?     , ?          , ?           , ?     , ?           ) "
                                 , (Isolate, Id_strain, Year, Month, Day, Host, Source, Institut   , Country_cod3, Region, Region_full ))
    Id_isolate=c.lastrowid
    c.execute("INSERT INTO isolate_seq (Id_isolate   , Id_seq) VALUES (?,?) "
                                     , (Id_isolate   , Id_seq))
    # Id_isolate_seq= c.fetchone()
    # Id_isolate_seq = Id_isolate_seq[0] if Id_isolate_seq else Id_isolate_seq

    if Id_algseq is None:
        #success = abnormal_row(c, row)
        c.execute("INSERT INTO pending_seq (Id_taxa, Name,      Id_seq) VALUES (?,?,?)",
                  (Id_taxa, MEGA_name, Id_seq))

        print("Abnormal row !!!!! ", MEGA_name, subtype, Str_name,
          "-------> Taxa:{0}, Alseq:{1}, Seq:{2}".format(Id_taxa, Id_algseq, Id_seq))
        success = False
    else:
        c.execute("INSERT INTO classified_seq (Id_taxa, Id_algseq) VALUES (?,?) "
                                            , (Id_taxa, Id_algseq)                 )
        # if c.rowcount == 0: return abnormal_row(c, row)

    rfs = []
    if Lu_Li         : rfs.append( 'Lu' )
    if Reference     : rfs.append( 'VR' )
    if Reference=='R': rfs.append( 'ICVT')
    print (rfs)
    for rf in rfs:
       print ('Add scheme: ',Id_taxa ,     rf                ,     Id_seq,   MEGA_name    )
       c.execute("INSERT INTO ref_seq (Id_taxa,     Id_ref_schema     ,      Id_seq,    name        ) "
               "             VALUES (?, (SELECT Id_ref_schema FROM ref_schema WHERE schema=?), ?, ?) "
                                  , (Id_taxa ,     rf  ,                     Id_seq,   MEGA_name    ))


    db.commit()
    return success

def parse_HEV_xlsm(db, file_name=None):
    if not file_name:
        file_name = filedialog.askopenfilename(filetypes=(("Excel files", "*.xlsm"), ("All files", "*.*")),
                                               defaultextension='fas',
                                               title='Select HEV isolate subtyping deta')
        if not file_name:
            return

    # self.refSeq.clear()
    print(file_name)

    wb = openpyxl.load_workbook(file_name)
    print(wb.sheetnames)

    ws = wb['Seq-class']   # 2 ('Seq-class')
    first = True
    error = False
    col=dict()
    for r in ws.iter_rows() :
        if first:
            for c in range(25):    # make a dic with the first 25 columns headers
                col[r[c].value]=c
            first = False
        else:
            error |= parse_row(db,r,col)
    if error:
        print('There were errors during parsing the Excel file !!!!!!!!!!!!!!')

def clean_parsed_Excel(db):
    c = db.cursor()
    c.execute("DELETE FROM strain")   # ??
    c.execute("DELETE FROM isolate")   # ??
    c.execute("DELETE FROM isolate_seq")   # ??
    c.execute("DELETE FROM pending_seq")   # ??
    c.execute("DELETE FROM classified_seq")   # ??
    db.commit()

def parseGB(db, GB_flat_file=None):
    '''
    Load and parse a GenBank sequence flat file
    :return:
    '''
    if not GB_flat_file:
        GB_flat_file = filedialog.askopenfilename(filetypes=(("Seq flat GB", "*.gb"), ("All files", "*.*") ),
                                                      title='Parse the GenBank sequences in flat format')
    print(GB_flat_file)

    c = db.cursor()

    # no real need to save this info
    c.execute("INSERT INTO seq_file (path, format) VALUES (?, 'GB_flat')", (GB_flat_file,))
    Id_file = c.lastrowid
    with open(GB_flat_file) as GB_flat:
      for record in GenBank.parse(GB_flat):
        Name     =  str(record.locus)   # todo: reparse
        strain   = ''
        isolate  = ''
        host     = ''
        country  = ''
        region   = ''
        collection_date = None
        source   = ''
        genotype = ''
        subtype  = None
        NCBI_TaxID = None

        #  http://biopython.org/DIST/docs/api/Bio.GenBank.Record-pysrc.html#Feature
        for feature in record.features:
            if feature.key == 'source':
                # http://biopython.org/DIST/docs/api/Bio.GenBank.Record.Qualifier-class.html
                for q in feature.qualifiers:

                    if q.key                                    == '/strain=':
                        strain = q.value[1:-1].strip()

                    elif q.key                                  == '/isolate=':
                        isolate = q.value[1:-1].strip()

                    elif q.key                                  == '/country=':
                        country = q.value[1:-1].split(':')
                        if len(country) > 1:
                            region = country[1].strip()  # ok?
                        country = country[0].strip()

                    elif q.key                                  == '/collection_date=':
                        collection_date = q.value[1:-1].strip()

                    elif q.key == '/source=' or q.key           == '/isolation_source=':
                        source = q.value[1:-1].strip()

                    elif q.key                                  == '/host=':
                        host = q.value[1:-1]

                    elif q.key                                  == '/db_xref=':      # example value= taxon:509628
                         val = q.value[1:-1]
                         m = val.split(':')
                         if m[0].startswith                         ('taxon'):
                             NCBI_TaxID = m[1].strip()

                    elif q.key                                  == '/note=':
                        val = q.value[1:-1]
                        print('Note=', val)
                        for n in val.split(';'):
                            m = n.lower().split()  #
                            if len(m) ==1 : m = m[0].split(':')
                            if m[0].startswith                      ('genotype'):
                                genotype = m[1].strip()
                            elif m[0].startswith                    ('subtype'):
                                subtype = m[1].strip()

        des = record.definition
        if 'isolate' in des:
            iso = des.split('isolate')[1].split(',')[0]
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
            st = des.split('strain')[1].split(',')[0]
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


        c.execute("INSERT INTO seq (Name,               Seq,     Len  ) "
                  "         VALUES (?,                  ?,       ?    )",
                                   (Name, str(record.sequence), len(record.sequence))    )
        Id_seq = c.lastrowid

        taxa= subtype if subtype else genotype
        if taxa:
            c.execute("SELECT Id_taxa FROM taxa WHERE taxa.Name=?", (taxa,))  # ?? Name UNIQUE ??
        else:
            c.execute("SELECT Id_taxa FROM taxa WHERE taxa.NCBI_TaxID=?", (NCBI_TaxID,)) # (record.organism,))
        Id_taxa = c.fetchone()
        Id_taxa = Id_taxa[0] if Id_taxa else Id_taxa

        # todo: parse location. Is unique?
        if not strain: strain = isolate
        if not isolate: isolate = strain
        c.execute("SELECT Id_strain FROM strain WHERE Name=?", (strain,))
        Id_strain = c.fetchone()
        if Id_strain:
            print('Existing Strain:', strain, Id_strain )
            c.execute("INSERT INTO strain (Name) VALUES (?) ", (strain,))
            Id_strain = c.lastrowid
            print('     new Id_strain:', Id_strain, )
        else:
            # print('New Strain:', Str_name, Id_strain)
            c.execute("INSERT INTO strain (Name) VALUES (?) ", (strain,))
            Id_strain = c.lastrowid
            # print('New Strain ID:', Id_strain)

        c.execute("SELECT iso3 FROM countries WHERE name=?", (country,))
        Country_cod3 = c.fetchone()
        Country_cod3 = Country_cod3[0] if Country_cod3 else None
        c.execute(
            "INSERT INTO isolate (Name   , Id_strain, year,            month, day, host, source, institution, country_iso3, region, region_full ) "
            "             VALUES (?      , ?        , ?   , ?    , ?  , ?   , ?     , ?          , ?           , ?     , ?           ) "
            ,                    (isolate, Id_strain, collection_date, None, None, host, source, None       , Country_cod3, None  , region      ))
        Id_isolate = c.lastrowid
        c.execute("INSERT INTO isolate_seq (Id_isolate   , Id_seq) VALUES (?,?) "
                  , (Id_isolate, Id_seq))
        # Id_isolate_seq= c.fetchone()
        # Id_isolate_seq = Id_isolate_seq[0] if Id_isolate_seq else Id_isolate_seq

        c.execute("INSERT INTO pending_seq (Name, Id_taxa, description, Id_isolate, Id_seq) VALUES (?,?,?,?,?)",
                  (Name, Id_taxa, record.definition, Id_isolate, Id_seq))

        db.commit()



if __name__ == '__main__':

    # exit(0)
    # """

    newly=True   # False
    country=True # False

    print('Creating db...')
    sdb = create(newly)

    if newly or country:
        read_country_codes(sdb)

    ref_name = "M73218"
    if newly:
        print('Parsing GB_flat_file...')
        parseGB(sdb, r'C:/Prog/HEV/data/temp/HEV-g3marked_all_2017-09-27.  577 seq.sequence.gb')
        print('Parsing the big alignment...')
        ID_align, ref = parse_full_fasta_Align(sdb, ref_name,'C:/Prog/HEV/alignment/HEV.fas')
        print(ref)
    else:
        print('Cleaning parsed Excel...')
        clean_parsed_Excel(sdb)
        ID_align=1      # ??

    print('Calcule reference positions...')
    ref = ref_pos(sdb, ID_align, ref_name) # , ref_name
    print(ref)

    print('Parse the Excel file table...')
    parse_HEV_xlsm(sdb, 'C:/Prog/HEV/data/temp/HEVsubtypingMEGAut - Kopie.xlsm')

    print('Done !!!')
    sdb.close()

    # """