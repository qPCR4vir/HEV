from tkinter import filedialog
import sqlite3
from Bio import SeqIO

# import tkinter

def read_create(sdb):
    with open("create_seq.sql") as dbcreate:
        sql_create = dbcreate.read()
    c = sdb.cursor()
    c.executescript(sql_create)


def create() -> sqlite3.Connection:
    sdb = sqlite3.connect("../data/temp/seq.db")
    read_create(sdb)
    return sdb


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
    ref_name = "M73218"
    ID_align, ref = parse_full_fasta_Align(sdb, ref_name)
    print(ref)

    ref = ref_pos(sdb, ID_align, ref_name) # , ref_name
    print(ref)
    sdb.close()
