    -- http://biosql.org/wiki/Schema_Overview
	-- https://github.com/biosql/biosql/blob/master/sql/biosqldb-sqlite.sql
	
	CREATE TABLE IF NOT EXISTS  seq_file
                     (
                       Id_file   INTEGER PRIMARY KEY AUTOINCREMENT,
                       path      TEXT    ,
                       format    TEXT
                     );

    -- seq    -- Any "experimental" sequence, for example: most GenBank entries (the sequence part).
    --           but also sequences obtained in the lab.
    CREATE TABLE IF NOT EXISTS  seq
                     (
                       Id_seq   INTEGER PRIMARY KEY AUTOINCREMENT,
                       Name     TEXT    UNIQUE,
                       Seq      TEXT,                  -- The original, preferable full sequence
                       -- Id_gen   INTEGER /* NOT NULL */ REFERENCES genomic_region,   --  (CG, ORF1, etc.)
                       Len      INT
                     );

    -- seq_file_pos
    --           Original file source position of the referenced sequence
    CREATE TABLE IF NOT EXISTS  seq_file_pos
                     (
                       Id_filep INTEGER PRIMARY KEY AUTOINCREMENT,
                       Id_seq   INTEGER REFERENCES seq,
                       Id_file  INTEGER REFERENCES seq_file,
                       file_beg INT,
                       file_end INT,      -- optional
                       time     INT       -- format ??
                     );

    -- bio_seq  -- A (preferably) complete sequence from a gen, organism, etc. For example the HEV complete genomes.
    -- temp: just a seq with Id_gen from a genomic_region with Name "CG".
    CREATE TABLE IF NOT EXISTS  bio_seq
           (
             Id_bio_seq       INTEGER PRIMARY KEY AUTOINCREMENT,
             Id_taxa          TEXT    UNIQUE,
             description      TEXT,
             locus            TEXT,
             Id_seq           INTEGER   REFERENCES seq
           );

    -- GB_seq  -- a simplified view of a GenBank entry sequence: todo use BioSQL
    CREATE TABLE IF NOT EXISTS  GB_seq
           (
             Id_GB_seq        INTEGER PRIMARY KEY AUTOINCREMENT,
             Acc              TEXT    UNIQUE,
             description      TEXT,
             locus            TEXT,
             Id_seq           INTEGER   REFERENCES seq
           );

    -- seq_fragment  -- an arbitrary fragment from some seq, for example: a BLAST  hit for which the
    --            original experimental sequence is known, or simple a region selected for phylogenetic analysis
    CREATE TABLE IF NOT EXISTS  seq_fragment
                     (
                       Id_frag INTEGER PRIMARY KEY NOT NULL,
                       Id_seq  INTEGER REFERENCES seq,
                       beg     INT,             -- relative to the original
                       end     INT
                     );
          
    -- partial_seq  -- Usually a temporal information construct containing the seq text,
    --                 for example from: a BLAST hit, or a region selected for phylogenetic analysis
    --       if the corresponding experimental seq is found it can be converted into a seq_fragment
    CREATE TABLE IF NOT EXISTS  partial_seq
                     (
                       Id_part INTEGER PRIMARY KEY NOT NULL,
                       Id_seq  INTEGER REFERENCES seq,  -- if known
                       beg     INT,             -- relative to the original
                       end     INT,
                       seq     TEXT              -- The aligned sequence, with internal gaps -, but without external
                     );
          
    -- align -- for axample a BLAST result from one query, or a "normal" Multialignment.
    --          if all the sequences are experimental (not necessary bio_seq) it is a full_align (don't BLAST!)
    CREATE TABLE IF NOT EXISTS    align
                     (
                       Id_align  INTEGER PRIMARY KEY AUTOINCREMENT,
                       Name      TEXT,
                       Id_file   INTEGER          REFERENCES seq_file,
                       Al_len    integer ,
                       Ref       text        -- if not NULL, the segested reference sequence
                     );
          
    -- aligned_seq
    CREATE TABLE IF NOT EXISTS  aligned_seq
                     (
                       Id_algseq  INTEGER PRIMARY KEY NOT NULL,
                       Id_align   INTEGER NOT NULL REFERENCES align,
                       Id_part    INTEGER NOT NULL REFERENCES seq, -- partial_seq,
                       Seq        TEXT,            -- The aligned sequence, with internal gaps
                       beg        INT,             -- relative to the alignment
                       end        INT
                     );
          
    -- strain
    CREATE TABLE IF NOT EXISTS   strain
           (
             Id_strain  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name       TEXT    UNIQUE
           );

    -- isolate
    CREATE TABLE IF NOT EXISTS   isolate
           (
             Id_isolate  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name        TEXT,
             Id_strain   INTEGER NOT NULL REFERENCES strain,
             idate        TEXT,
             year        INT,
             month       INT,
             day         INT,
             Id_host     INTEGER NOT NULL,     -- original taxa
             Id_source   INTEGER NOT NULL,
             Id_author   INTEGER NOT NULL,
             Id_institution INTEGER NOT NULL,
             Id_location    INTEGER NOT NULL,
             coordinate  TEXT
           );

    -- author
    CREATE TABLE IF NOT EXISTS  author
           (
             Id_author  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name       TEXT    ,
             middle     TEXT ,
             surname    TEXT ,
             title      TEXT,
             genre      TEXT,
             Id_contact INTEGER NOT NULL
           );

    -- institution
    CREATE TABLE IF NOT EXISTS  institution
           (
             Id_institution  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name            TEXT    ,
             category        TEXT ,                    -- Institute, Department, University, Laboratory, Center, etc.
             parent          INTEGER NOT NULL,
             Id_contact      INTEGER NOT NULL
           );

    -- affiliation
    CREATE TABLE IF NOT EXISTS  affiliation
           (
             Id_author        INTEGER NOT NULL,
             Id_institution   INTEGER NOT NULL,
             beg_year      INT,
             beg_month     INT,
             end_year      INT,
             end_month     INT,
             Id_contact    INTEGER NOT NULL,
             position      TEXT
           );

    -- genomic_region
    CREATE TABLE IF NOT EXISTS  genomic_region
           (
             Id_gen   INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    ,
             Id_taxa  INTEGER NOT NULL REFERENCES taxa,
             beg      INT,
             end      INT,
             Id_seq   INTEGER   REFERENCES seq   -- just to map the positions, preferably a bio_seq
           );

    -- source
    CREATE TABLE IF NOT EXISTS   source
           (
             Id_host        INTEGER NOT NULL,     -- original taxa
             artificial     TEXT ,                -- NULL, clone, cell, patent, etc.
             organ          TEXT,                 -- lever, serum, faces, etc.
             environment    TEXT                  -- NULL, seawater,  etc.
            );

    -- contact_adress
    CREATE TABLE IF NOT EXISTS  contact_adress
           (
             e_mail    TEXT,
             URL       TEXT,
             mail      TEXT,
             telephone TEXT,
             country   INTEGER NOT NULL,
             FAX       TEXT
           );

    -- country
    CREATE TABLE IF NOT EXISTS  country
           (
             Id_country INTEGER PRIMARY KEY AUTOINCREMENT,
             Name       TEXT  UNIQUE  ,
             idx        TEXT  UNIQUE ,
             region     TEXT,
             continent  TEXT
           );

    -- taxa
    CREATE TABLE IF NOT EXISTS  taxa
           (
             Id_taxa     INTEGER PRIMARY KEY AUTOINCREMENT,
             Name        TEXT    ,
             vulgar      TEXT ,
             parent      INTEGER NOT NULL,
             Id_category INTEGER NOT NULL,
             NCBI        TEXT
           );

    -- taxa_category
    CREATE TABLE IF NOT EXISTS  taxa_category
           (
             Id_category   INTEGER PRIMARY KEY AUTOINCREMENT,
             Name          TEXT    UNIQUE,
             parent        INTEGER NOT NULL,
             NCBI          TEXT
           );

    -- to_excel  VIEW
    CREATE VIEW  files AS SELECT  path, format FROM seq_file;
    CREATE VIEW  to_excel AS SELECT  Name, Len FROM seq;
    CREATE VIEW  all_frag AS SELECT  Name, Len, beg, end FROM seq, aligned_seq ON seq.Id_seq = aligned_seq.Id_part ;
            /*
             Seq_Name  ,
             parent ,
             NCBI
           )*/

