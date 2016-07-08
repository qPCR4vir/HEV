import tkinter
from tkinter import filedialog
import sqlite3
sdb = None



def create() -> sqlite3.Connection:
    global sdb
    sdb = sqlite3.connect("../seq.db")
    c = sdb.cursor()

    # seq_file
    c.execute("""CREATE TABLE   seq_file
                     (
                       Id     INTEGER PRIMARY KEY AUTOINCREMENT,
                       path   TEXT    ,
                       format TEXT,
                     )
              """)

    # seq
    c.execute("""CREATE TABLE   seq
                     (
                       Id    INTEGER PRIMARY KEY AUTOINCREMENT,
                       Name  TEXT    UNIQUE,
                       Seq   TEXT,                  -- The original, preferable full sequence
                       gen   INTEGER NOT NULL,      -- Genomic region (CG, ORF1, etc.)
                       Len   INT,
                       file  INTEGER,
                       file_beg INT,
                       file_end INT,
                       time  INT
                     )
              """)

    # seq_fragment
    c.execute("""CREATE TABLE   seq_fragment
                     (
                       Id      INTEGER PRIMARY KEY NOT NULL,
                       seq_ori INTEGER  ,
                       beg     INT,             -- relative to the original
                       end     INT,
                       seq     TXT              -- The aligned sequence, with internal gaps -, but without external
                     )
              """)

    # align
    c.execute(""""CREATE TABLE    align
                     (
                       Id    INTEGER PRIMARY KEY AUTOINCREMENT,
                       Name  TEXT
                     )
              """)

    # aligned_seq
    c.execute("""CREATE TABLE   aligned_seq
                     (
                       Id      INTEGER PRIMARY KEY NOT NULL,
                       seq_ori INTEGER ,
                       seq_frg INTEGER ,
                       align   INTEGER NOT NULL,
                       beg     INT,             -- relative to the alignment
                       end     INT,
                       seq     TXT              -- The aligned sequence, with internal gaps -, but without external
                     )
              """)

    # strain
    c.execute("""CREATE TABLE   strain
           (
             Id    INTEGER PRIMARY KEY AUTOINCREMENT,
             Name  TEXT    UNIQUE
           )
    """)

    # isolate
    c.execute("""CREATE TABLE   isolate
           (
             Id     INTEGER PRIMARY KEY AUTOINCREMENT,
             Name   TEXT,
             strain INTEGER NOT NULL,
             date   TEXT,
             year   INT,
             month  INT,
             day    INT,
             host   INTEGER NOT NULL,     -- original taxa
             source INTEGER NOT NULL,
             author INTEGER NOT NULL,
             institution INTEGER NOT NULL,
             location    INTEGER NOT NULL,
             coordinate  TXT,
           )
    """)

    # author
    c.execute("""CREATE TABLE   author
           (
             Id       INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    ,
             middle   TEXT ,
             surname  TEXT ,
             title    TEXT,
             genre    TEXT,
             contact  INTEGER NOT NULL
           )
    """)

    # institution
    c.execute("""CREATE TABLE   institution
           (
             Id       INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    ,
             category TEXT ,
             parent   INTEGER NOT NULL,
             contact  INTEGER NOT NULL,
           )
    """)

    # affiliation
    c.execute("""CREATE TABLE   affiliation
           (
             author        INTEGER NOT NULL,
             institution   INTEGER NOT NULL,
             beg_year      INT,
             beg_month     INT,
             end_year      INT,
             end_month     INT,
             contact       INTEGER NOT NULL,
             position      TEXT

           )
    """)

    # genomic_region
    c.execute("""CREATE TABLE   genomic_region
           (
             Id    INTEGER PRIMARY KEY AUTOINCREMENT,
             Name  TEXT    ,
             taxa  INTEGER NOT NULL,
             beg   INT,
             end   INT
           )
    """)

    # source
    c.execute("""CREATE TABLE   source
           (
             host        INTEGER NOT NULL,     -- original taxa
             artificial  TEXT ,                -- NULL, clone, cell, patent, etc.
             organ       TEXT,                 -- lever, serum, faces, etc.
             environment TEXT                  -- NULL, seawater,  etc.
            )
    """)

    # contact_adress
    c.execute("""CREATE TABLE   contact_adress
           (
             e_mail    TEXT,
             URL       TEXT,
             mail      TEXT,
             telephone TEXT,
             country   INTEGER NOT NULL,
             FAX       TEXT
           )
    """)

    # country
    c.execute("""CREATE TABLE   country
           (
             Id        INTEGER PRIMARY KEY AUTOINCREMENT,
             Name      TEXT  UNIQUE  ,
             idx       TEXT  UNIQUE ,
             region    TEXT,
             continent TEXT
           )
    """)

    # taxa
    c.execute("""CREATE TABLE   taxa
           (
             Id       INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    ,
             vulgar   TEXT ,
             parent   INTEGER NOT NULL,
             category INTEGER NOT NULL,
             NCBI     TXT
           )
    """)

    # taxa_category
    c.execute("""CREATE TABLE   taxa_category
           (
             Id       INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    UNIQUE,
             parent   INTEGER NOT NULL,
             NCBI     TXT
           )
    """)

    # big VIEW
    c.execute("""CREATE VIEW TABLE   big   AS SELECT
           (
             Id       INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    UNIQUE,
             parent   INTEGER NOT NULL,
             NCBI     TXT
           )
    """)

    return sdb


def parseAlign( file_name=None):
    if not file_name:
        file_name = filedialog.askopenfilename(filetypes=(("fasta aligment", "*.fas"), ("All files", "*.*")),
                                               defaultextension='fas',
                                               title='Select Master Alignment')
        if not file_name:
            return

    # self.refSeq.clear()
    print(file_name)
    c = sdb.cursor()
    c.execute("INSERT INTO align VALUES (?,?)", (1, file_name))
    with open(file_name) as align_file:
        seq_name = ''
        seq_num = 0
        for line in align_file.readlines():
            # print(line)
            if line[0] == '>':
                # todo: end previous seq !?
                seq_num += 1
                seq_name = line[1:].rstrip()
            else:
                ln = len(line)
                #self.refLen = max(ln, self.refLen)
                seq_beg = 0
                seq_end = ln - 2
                while seq_beg < ln:
                    if line[seq_beg] == '-':
                        seq_beg += 1
                    else:
                        break  # todo :  check it is a valid base not line end???
                while seq_end > seq_beg:
                    if line[seq_end] == '-':
                        seq_end -= 1
                    else:
                        break  # todo :  check it is a valid base not line end???
                seq = line[seq_beg: seq_end+1]
                print('Seq: ' + seq_name + str((seq_beg, seq_end)) + ": " + seq)

                # CREATE TABLE seq(                            Id INT, Name TEXT, Seq TEXT, Len INT)")
                c.execute("INSERT INTO seq VALUES (?,?,?,?)", (seq_num, seq_name, seq, len(seq)))

                # self.refSeq[seq_name] = Seq_pos(seq_beg, seq_end)
    sdb.commit()
    #self.ID_original.clear()
    #self.ID_original.add('\n'.join(self.refSeq.keys()))


if __name__=='__main__':
    create()
    parseAlign()
    sdb.close()