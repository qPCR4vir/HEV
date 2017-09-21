from tkinter import filedialog
import sqlite3
from Bio import SeqIO

# import tkinter

def create() -> sqlite3.Connection:
    sdb = sqlite3.connect("../data/temp/seq.db")
    read_create(sdb)
    return sdb


def read_create(sdb):
    with open("create_seq.sql") as dbcreate:
        sql_create = dbcreate.read()
    c = sdb.cursor()
    c.executescript(sql_create)


def add_def_taxa(sdb):
    c = sdb.cursor()
    #  ----------------------------------------------------------------------------------------------------
    c.execute("INSERT INTO taxa_rank (Name,           NCBI          )"     
              "               VALUES ('superkingdom', 'superkingdom')"   )
    parent_rank = c.lastrowid

    c.execute("INSERT INTO taxa      (Name     ,  vulgar  , Id_rank   , NCBI_TaxID   )"
              "               VALUES ('Viridae', 'viruses', ?         , '10239'      )",
                                     (                     parent_rank,              ) )
    parent_taxa = c.lastrowid

    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa_rank (Name     , parent    , NCBI          )"
              "               VALUES ('no rank', ?         ,'no rank'      )"   ,
                                     (          parent_rank,               )    )
    parent_rank = c.lastrowid

    c.execute("INSERT INTO taxa      (Name           ,  vulgar  , Id_rank    , parent,  NCBI_TaxID   )"
              "               VALUES ('ssRNA viruses', 'viruses', ?          , ?     , '439488'      )",
                                     (                            parent_rank, parent_taxa           ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa_rank (Name     , parent    , NCBI          )"
              "               VALUES (NULL, ?         ,'no rank'      )"   ,
                                     (          parent_rank,               )    )
    parent_rank = c.lastrowid

    c.execute("INSERT INTO taxa      (Name                                         ,  vulgar  , Id_rank, parent , NCBI_TaxID )"
              "               VALUES ('ssRNA positive-strand viruses, no DNA stage', 'viruses', ?       , ?     , '35278'    )",
                                     (                                                       parent_rank, parent_taxa            ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa_rank (Name     , parent    , NCBI          )"
              "               VALUES ('family', ?          ,'family'       )"   ,
                                     (          parent_rank,               )    )
    parent_rank = c.lastrowid

    c.execute("INSERT INTO taxa      (Name         ,  vulgar , Id_rank    , parent,  NCBI_TaxID   )"
              "               VALUES ('Hepeviridae', 'HEV'   , ?          , ?     , '291484'      )",
                                     (                        parent_rank , parent_taxa           ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa_rank (Name     , parent    , NCBI          )"
              "               VALUES ('genus'  , ?         ,'genus'        )"   ,
                                     (          parent_rank,               )    )
    parent_rank = c.lastrowid

    c.execute("INSERT INTO taxa      (Name            ,  vulgar          , Id_rank     , parent ,  NCBI_TaxID   )"
              "               VALUES ('Orthohepevirus', 'Orthohepevirus' , ?           , ?      , '1678141'     )",
                                     (                                     parent_rank , parent_taxa            ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa_rank (Name     , parent    , NCBI          )"
              "               VALUES ('species', ?         ,'species'      )"   ,
                                     (          parent_rank,               )    )
    parent_rank_s = c.lastrowid

    c.execute("INSERT INTO taxa      (Name              ,  vulgar            , Id_rank     , parent ,  NCBI_TaxID   )"
              "               VALUES ('Orthohepevirus A', 'Orthohepevirus A' , ?           , ?      , '1678143'     )",
                                     (                                         parent_rank_s , parent_taxa            ) )
    parent_taxa_A = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa_rank (Name      , parent    , NCBI          )"
              "               VALUES ('genotype', ?         ,'no rank'      )"   ,
                                     (          parent_rank_s,               )    )
    parent_rank_g = c.lastrowid

    c.execute("INSERT INTO taxa      (Name ,  vulgar  , Id_rank     , parent ,  NCBI_TaxID   )"
              "               VALUES ('1'  , 'HEV-g1' , ?           , ?      , '185579'     )",
                                     (                 parent_rank_g , parent_taxa_A         ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar  , Id_rank     , parent ,  NCBI_TaxID   )"
              "               VALUES ('2'  , 'HEV-g2' , ?           , ?      , ''     )",
                                     (                 parent_rank_g , parent_taxa_A         ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar  , Id_rank     , parent ,  NCBI_TaxID   )"
              "               VALUES ('3'  , 'HEV-g3' , ?           , ?      , '509628'     )",
                                     (                 parent_rank_g , parent_taxa_A         ) )
    taxa_g3 = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar  , Id_rank     , parent ,  NCBI_TaxID   )"
              "               VALUES ('4'  , 'HEV-g4' , ?           , ?      , '185580'     )",
                                     (                 parent_rank_g , parent_taxa_A         ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar  , Id_rank     , parent ,  NCBI_TaxID   )"
              "               VALUES ('5'  , 'HEV-g5' , ?           , ?      , ''     )",
                                     (                 parent_rank_g , parent_taxa_A         ) )
    parent_taxa = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa_rank (Name      , parent    , NCBI          )"
              "               VALUES ('subtype' , ?         ,'no rank'      )"   ,
                                     (          parent_rank_g,               )    )
    rank_subtype = c.lastrowid


    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3a'  , 'HEV-g3a' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3a = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3b'  , 'HEV-g3b' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3b = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3c'  , 'HEV-g3c' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3c = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3d'  , 'HEV-g3d' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3d = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3e'  , 'HEV-g3e' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3e = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3f'  , 'HEV-g3f' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3f = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3g'  , 'HEV-g3g' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3g = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3h'  , 'HEV-g3h' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3h = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3i'  , 'HEV-g3i' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3i = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3j'  , 'HEV-g3j' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3j = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------

    c.execute("INSERT INTO taxa      (Name ,  vulgar    , Id_rank      , parent ,  NCBI_TaxID   )"
              "               VALUES ('3k'  , 'HEV-g3k' , ?            , ?      , ''     )",
                                     (                    rank_subtype , taxa_g3         ) )
    taxa_g3k = c.lastrowid
    #  ----------------------------------------------------------------------------------------------------



def parse_full_fasta_Align( sdb, ref_seq = None, file_name=None):
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

    c = sdb.cursor()
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
