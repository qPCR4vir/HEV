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
                       Id_file   INTEGER PRIMARY KEY AUTOINCREMENT,
                       path      TEXT    ,
                       format    TEXT,
                     );
              """)

    # seq    -- Any "experimental" sequence, for example: most GenBank entries (the sequence part).
    #           but also sequences obtained in the lab.
    c.execute("""CREATE TABLE   seq
                     (
                       Id_seq   INTEGER PRIMARY KEY AUTOINCREMENT,
                       Name     TEXT    UNIQUE,
                       Seq      TEXT,                  -- The original, preferable full sequence
                       Id_gen   INTEGER NOT NULL REFERENCES genomic_region,   --  (CG, ORF1, etc.)
                       Len      INT,
                       Id_file  INTEGER REFERENCES seq_file,
                       file_beg INT,
                       file_end INT,
                       time     INT
                     );
              """)

    # bio_seq  -- A complete sequence from a gen, organism, etc. For example the HEV complete genomes.
    # temp: just a seq with Id_gen from a genomic_region with Name "CG".
    c.execute("""CREATE TABLE   bio_seq
           (
             Id_bio_seq       INTEGER PRIMARY KEY AUTOINCREMENT,
             Id_taxa          TEXT    UNIQUE,
             description      TEXT,
             locus            TEXT,
             Id_seq           INTEGER   REFERENCES seq
           );
    """)

    # GB_seq  -- a simplified view of a GenBank entry sequence: todo use BioSQL
    c.execute("""CREATE TABLE   GB_seq
           (
             Id_GB_seq        INTEGER PRIMARY KEY AUTOINCREMENT,
             Acc              TEXT    UNIQUE,
             description      TEXT,
             locus            TEXT,
             Id_seq           INTEGER   REFERENCES seq
           );
    """)

    # seq_fragment  -- just an arbitrary fragment from some seq, for example: a BLAST  hit,
    #                  or a region selected for phylogenetic analysis
    c.execute("""CREATE TABLE   seq_fragment
                     (
                       Id_frag INTEGER PRIMARY KEY NOT NULL,
                       Id_seq  INTEGER REFERENCES seq,
                       beg     INT,             -- relative to the original
                       end     INT,
                       seq     TXT              -- The aligned sequence, with internal gaps -, but without external
                     );
              """)

    # align -- for axample a BLAST result from one query, or a "normal" Multialignment.
    #          if all the sequences are experimental (not necessary bio_seq) it is a full_align (don't BLAST!)
    c.execute(""""CREATE TABLE    align
                     (
                       Id_align  INTEGER PRIMARY KEY AUTOINCREMENT,
                       Name      TEXT
                       Id_file   INTEGER          REFERENCES seq_file,
                     );
              """)

    # aligned_seq
    c.execute("""CREATE TABLE   aligned_seq
                     (
                       Id_algseq  INTEGER PRIMARY KEY NOT NULL,
                       Id_seq     INTEGER          REFERENCES seq,
                       Id_frag    INTEGER          REFERENCES seq_fragment,
                       Id_align   INTEGER NOT NULL REFERENCES align,
                       beg        INT,             -- relative to the alignment
                       end        INT,
                       seq        TXT              -- The aligned sequence, with internal gaps -, but without external
                     )
              """)

    # strain
    c.execute("""CREATE TABLE   strain
           (
             Id_strain  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name       TEXT    UNIQUE
           )
    """)

    # isolate
    c.execute("""CREATE TABLE   isolate
           (
             Id_isolate  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name        TEXT,
             Id_strain   INTEGER NOT NULL,
             date        TEXT,
             year        INT,
             month       INT,
             day         INT,
             Id_host     INTEGER NOT NULL,     -- original taxa
             Id_source   INTEGER NOT NULL,
             Id_author   INTEGER NOT NULL,
             Id_institution INTEGER NOT NULL,
             Id_location    INTEGER NOT NULL,
             coordinate  TXT,
           )
    """)

    # author
    c.execute("""CREATE TABLE   author
           (
             Id_author  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name       TEXT    ,
             middle     TEXT ,
             surname    TEXT ,
             title      TEXT,
             genre      TEXT,
             Id_contact INTEGER NOT NULL
           )
    """)

    # institution
    c.execute("""CREATE TABLE   institution
           (
             Id_institution  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name            TEXT    ,
             category        TEXT ,                    -- Institute, Department, University, Laboratory, Center, etc.
             parent          INTEGER NOT NULL,
             Id_contact      INTEGER NOT NULL,
           )
    """)

    # affiliation
    c.execute("""CREATE TABLE   affiliation
           (
             Id_author        INTEGER NOT NULL,
             Id_institution   INTEGER NOT NULL,
             beg_year      INT,
             beg_month     INT,
             end_year      INT,
             end_month     INT,
             Id_contact    INTEGER NOT NULL,
             position      TEXT

           )
    """)

    # genomic_region
    c.execute("""CREATE TABLE   genomic_region
           (
             Id_gen   INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    ,
             Id_taxa  INTEGER NOT NULL REFERENCES taxa,
             beg      INT,
             end      INT,
             Id_seq   INTEGER   REFERENCES seq,   -- just to map the positions, preferably a bio_seq
           )
    """)

    # source
    c.execute("""CREATE TABLE   source
           (
             Id_host        INTEGER NOT NULL,     -- original taxa
             artificial     TEXT ,                -- NULL, clone, cell, patent, etc.
             organ          TEXT,                 -- lever, serum, faces, etc.
             environment    TEXT                  -- NULL, seawater,  etc.
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
             Id_country INTEGER PRIMARY KEY AUTOINCREMENT,
             Name       TEXT  UNIQUE  ,
             idx        TEXT  UNIQUE ,
             region     TEXT,
             continent  TEXT
           )
    """)

    # taxa
    c.execute("""CREATE TABLE   taxa
           (
             Id_taxa     INTEGER PRIMARY KEY AUTOINCREMENT,
             Name        TEXT    ,
             vulgar      TEXT ,
             parent      INTEGER NOT NULL,
             Id_category INTEGER NOT NULL,
             NCBI        TXT
           )
    """)

    # taxa_category
    c.execute("""CREATE TABLE   taxa_category
           (
             Id_category   INTEGER PRIMARY KEY AUTOINCREMENT,
             Name          TEXT    UNIQUE,
             parent        INTEGER NOT NULL,
             NCBI          TXT
           )
    """)

    # to_excel  VIEW
    c.execute("""CREATE VIEW TABLE   to_excel   AS SELECT
           (
             -- Id       INTEGER PRIMARY KEY AUTOINCREMENT,
             Seq_Name     TEXT    UNIQUE,
             parent   INTEGER NOT NULL,
             NCBI     TXT
           )
    """)

    # to_doc


    return sdb


def parse_full_Align( file_name=None):
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
    c.execute("INSERT INTO seq_file (path) VALUES (?)", (file_name,))
    Id_file = c.lastrowid()
    c.execute("INSERT INTO align (Id_file, Name) VALUES (?, ?)", (Id_file, file_name))
    Id_align = c.lastrowid()

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
                c.execute("INSERT INTO seq (Id, Name, Seq, Len) VALUES (?,?,?,?)", (seq_num, seq_name, seq, len(seq)))

                # self.refSeq[seq_name] = Seq_pos(seq_beg, seq_end)
    sdb.commit()
    #self.ID_original.clear()
    #self.ID_original.add('\n'.join(self.refSeq.keys()))


if __name__=='__main__':
    create()
    parseAlign()
    sdb.close()