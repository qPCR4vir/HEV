from tkinter import filedialog
import sqlite3
from Bio import SeqIO

# import tkinter


class CreateTaxa:
    def __init__(self, db):
        self.db=db
        self.c = sdb.cursor()

    def root_rank (self, name, NCBI = None):        #  todo: NCBI ? add kingdom ! determine type of rank.
        self.c.execute("INSERT INTO taxa_rank (Name, NCBI )"
                       "               VALUES (?   , ?    )",
                       (                       name, NCBI )   )
        self.r_rank = self.c.lastrowid
        return self.r_rank

    def root_taxa (self, Name ,  vulgar, NCBI_TaxID):
        self.c.execute("INSERT INTO taxa      (Name , vulgar, Id_rank    , NCBI_TaxID  )"
                       "               VALUES ( ?   , ?     , ?          , ?           )",
                                              (Name , vulgar, self.r_rank, NCBI_TaxID  ) )
        self.r_taxa = self.c.lastrowid
        return self.r_taxa

    def rank (self, name, parent_rank  , NCBI = None):
        self.c.execute("INSERT INTO taxa_rank (Name, parent      , NCBI )"
                       "               VALUES (?   , ?           , ?    )",
                       (                       name, parent_rank , NCBI )   )
        return self.c.lastrowid

    def taxa (self, Name, vulgar, rank    , parent_taxa,  NCBI_TaxID  = None):
        self.c.execute("INSERT INTO taxa (Name,  vulgar, Id_rank, parent     , NCBI_TaxID )"
                       "          VALUES (?   , ?      , ?      , ?          , ?          )",
                       (                  Name, vulgar , rank   , parent_taxa, NCBI_TaxID )   )
        return self.c.lastrowid


def create() -> sqlite3.Connection:
    sdb = sqlite3.connect("../data/temp/seq.db")
    read_create(sdb)
    return sdb


def read_create(sdb):
    with open("create_seq.sql") as dbcreate:
        sql_create = dbcreate.read()
    c = sdb.cursor()
    c.executescript(sql_create)


def add_def_taxa(db):
    ct= CreateTaxa(db)

    ct.root_rank('superkingdom', 'superkingdom')
    ct.root_taxa('Viridae', 'viruses', '10239')

    r   = ct.rank('no rank', ct.r_rank, 'no rank')
    t   = ct.taxa('ssRNA viruses', 'viruses', r, ct.r_taxa, '439488')

    r   = ct.rank('', r, 'no rank')
    t   = ct.taxa('ssRNA positive-strand viruses, no DNA stage', 'viruses', r, t, '35278')

    r   = ct.rank('family', r, 'family')
    t   = ct.taxa('Hepeviridae', 'HEV'  , r, t, '291484')

    rG  = ct.rank('genus', r, 'genus')
    tO  = ct.taxa('Orthohepevirus', 'Orthohepevirus'  , rG, t, '1678141')

    rs  = ct.rank('species', rG, 'species')
    tA  = ct.taxa('Orthohepevirus A', 'Orthohepevirus A'  , rs, tO, '1678143')

    rg  = ct.rank('genotype', rs, 'genotype')
    g1  = ct.taxa('1', 'HEV-g1'  , rg, tA, '185579')
    g2  = ct.taxa('2', 'HEV-g2'  , rg, tA )
    g3  = ct.taxa('3', 'HEV-g3'  , rg, tA, '509628')
    g4  = ct.taxa('4', 'HEV-g4'  , rg, tA, '185580')
    g5  = ct.taxa('5', 'HEV-g5'  , rg, tA )
    g6  = ct.taxa('6', 'HEV-g6'  , rg, tA )
    g7  = ct.taxa('7', 'HEV-g7'  , rg, tA )

    rsubt = ct.rank('subtype', rg, 'no rank')
    
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

    g3a   = ct.taxa('3a', 'HEV-g3a'  , rsubt, g3)
    g3b   = ct.taxa('3b', 'HEV-g3b'  , rsubt, g3)
    g3c   = ct.taxa('3c', 'HEV-g3c'  , rsubt, g3)
    g3d   = ct.taxa('3d', 'HEV-g3d'  , rsubt, g3)
    g3e   = ct.taxa('3e', 'HEV-g3e'  , rsubt, g3)
    g3f   = ct.taxa('3f', 'HEV-g3f'  , rsubt, g3)
    g3g   = ct.taxa('3g', 'HEV-g3g'  , rsubt, g3)
    g3h   = ct.taxa('3h', 'HEV-g3h'  , rsubt, g3)
    g3i   = ct.taxa('3i', 'HEV-g3i'  , rsubt, g3)
    g3j   = ct.taxa('3j', 'HEV-g3j'  , rsubt, g3)
    g3k   = ct.taxa('3k', 'HEV-g3k'  , rsubt, g3)

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
    ref = None
    for seq_record in SeqIO.parse(file_name, "fasta"):
        # print(seq_record.id, len(seq_record) )
        if ref_seq == seq_record.id :
            sr = 0
            ref = []
            for b in seq_record.seq:
                if (b != '-'): sr += 1
                ref.append(sr)
        ln = len(seq_record.seq)
        if max_len < ln: max_len = ln
        seq_beg = 0
        seq_end = ln - 1
        while seq_beg < ln:
            if seq_record.seq[seq_beg] == '-':
                seq_beg += 1
            else:
                break  # todo :  check it is a valid base not line end???
        while seq_end > seq_beg:
            if seq_record.seq[seq_end] == '-':
                seq_end -= 1
            else:
                break  # todo :  check it is a valid base not line end???

        seq = seq_record.seq[seq_beg: seq_end + 1]

        exp_seq = ''.join([base for base in seq if base != '-'])

        c.execute("UPDATE align SET Al_len = ? WHERE Id_align=?", (max_len, Id_align) )

        c.execute("INSERT INTO seq (Name,               Seq,     Len  ) "
                  "         VALUES (?,                  ?,       ?    )",
                                   (str(seq_record.id), exp_seq, len(exp_seq))    )
        Id_part = c.lastrowid

        c.execute("INSERT INTO aligned_seq (Id_align, Id_part, Seq,      beg,     end  ) "
                  "                 VALUES (?,        ?,       ?,        ?,       ?    )",
                                           (Id_align, Id_part, str(seq), seq_beg, seq_end )    )

    sdb.commit()
    return Id_align , ref


def ref_pos(sdb, ID_align, seq_name=None):
    c = sdb.cursor()
    c.execute("SELECT Al_len, Ref FROM align WHERE Id_align=? ", (ID_align,  ))
    Al_len, Ref = c.fetchone()
    if not seq_name: seq_name = Ref

    c.execute("SELECT aligned_seq.Seq, beg, end FROM aligned_seq, Seq ON Id_part=Id_seq WHERE Id_align=? AND Name=?", (ID_align, seq_name))
    Seq, beg, end = c.fetchone()
    sr=0
    ref = [sr]*beg
    for b in Seq:
        if (b != '-'): sr+=1
        ref.append(sr)
    ref += [sr]*(Al_len - end-1)
    return ref



if __name__ == '__main__':
    sdb = create()
    add_def_taxa(sdb)
    ref_name = "M73218"
    ID_align, ref = parse_full_fasta_Align(sdb, ref_name)
    print(ref)

    ref = ref_pos(sdb, ID_align, ref_name) # , ref_name
    print(ref)
    sdb.close()
