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
                       Id_seq   INTEGER REFERENCES seq(Id_seq),
                       Id_file  INTEGER REFERENCES seq_file(Id_file),
                       file_beg INT,
                       file_end INT,      -- optional
                       time     INT       -- format ??
                     );

    -- bio_seq  -- A (preferably) complete sequence from a gen, organism, etc. For example the HEV complete genomes.
    -- temp: just a seq with Id_gen from a genomic_region with Name "CG".
    CREATE TABLE IF NOT EXISTS  bio_seq
           (
             Id_bio_seq       INTEGER PRIMARY KEY AUTOINCREMENT,
             Id_taxa          INTEGER   REFERENCES taxa(Id_taxa),
             description      TEXT,
             locus            TEXT,
             Id_seq           INTEGER   REFERENCES seq(Id_seq)
           );

    CREATE TABLE IF NOT EXISTS  isolate_seq
           (
             Id_isolate_seq   INTEGER   PRIMARY KEY AUTOINCREMENT,
             Id_isolate       INTEGER   REFERENCES isolate(Id_isolate),
             Id_seq           INTEGER   REFERENCES seq(Id_seq)
           );

    CREATE TABLE IF NOT EXISTS  classified_seq
           (
             Id_clas_seq      INTEGER PRIMARY KEY AUTOINCREMENT,
             Id_taxa          INTEGER   REFERENCES taxa,          -- the finest available classification
             -- description      TEXT,                               --  ??
             -- Id_isolate       INTEGER   REFERENCES isolate    ,
             Id_algseq        INTEGER  NOT NULL REFERENCES aligned_seq(Id_algseq)
           );
           );


    CREATE TABLE IF NOT EXISTS  pending_seq
           (
             Id_pend_seq      INTEGER   PRIMARY KEY AUTOINCREMENT,
             Name             TEXT,
             Id_taxa          INTEGER   REFERENCES taxa(Id_taxa),          -- the finest available classification
             description      TEXT,                               --  ??
             Id_isolate       INTEGER   REFERENCES isolate(Id_isolate)    ,
             Id_algseq        INTEGER   REFERENCES aligned_seq(Id_algseq),
             Id_seq           INTEGER   REFERENCES seq(Id_seq)
           );

   -- GB_seq  -- a simplified view of a GenBank entry sequence: todo use BioSQL
    CREATE TABLE IF NOT EXISTS  GB_seq
           (
             Id_GB_seq        INTEGER PRIMARY KEY AUTOINCREMENT,
             Acc              TEXT    UNIQUE,
             description      TEXT,
             locus            TEXT,
             Id_seq           INTEGER   REFERENCES seq(Id_seq)
           );

    -- seq_fragment  -- an arbitrary fragment from some seq, for example: a BLAST  hit for which the
    --            original experimental sequence is known, or simple a region selected for phylogenetic analysis
    CREATE TABLE IF NOT EXISTS  seq_fragment
                     (
                       Id_frag INTEGER PRIMARY KEY NOT NULL,
                       Id_seq  INTEGER REFERENCES seq(Id_seq),
                       pbeg     INT,             -- relative to the original
                       pend     INT
                     );
          
    -- partial_seq  -- Usually a temporal information construct containing the seq text,
    --                 for example from: a BLAST hit, or a region selected for phylogenetic analysis
    --       if the corresponding experimental seq is found it can be converted into a seq_fragment
    CREATE TABLE IF NOT EXISTS  partial_seq
                     (
                       Id_part INTEGER PRIMARY KEY NOT NULL,
                       Id_seq  INTEGER REFERENCES seq(Id_seq),  -- if known
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
                       Id_file   INTEGER          REFERENCES seq_file(Id_file),
                       Al_len    integer ,
                       Ref       text        -- if not NULL, the suggested reference sequence
                     );
          
    -- aligned_seq
    CREATE TABLE IF NOT EXISTS  aligned_seq
                     (
                       Id_algseq  INTEGER PRIMARY KEY NOT NULL,
                       Id_align   INTEGER NOT NULL REFERENCES align,
                       Id_part    INTEGER NOT NULL REFERENCES seq(Id_seq), -- partial_seq,
                       Seq        TEXT,            -- The aligned sequence, with internal gaps
                       pbeg        INT,             -- relative to the alignment
                       pend        INT
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
             Id_strain   INTEGER NOT NULL REFERENCES strain(Id_strain),
             -- idate       TEXT,                -- ?
             year        INT,
             month       INT,
             day         INT,
             host            TEXT,      -- todo: Id_host     INTEGER, --NOT NULL,     -- original taxa
             source          TEXT,      -- todo: Id_source   INTEGER, --NOT NULL,
             Id_author       INTEGER,   -- NOT NULL,
             institution     TEXT,      -- todo: Id_institution  INTEGER, --NOT NULL,
             Id_location     INTEGER REFERENCES location(Id_location), -- NOT NULL,
             Id_country_cod  TEXT REFERENCES countries(iso3)
           );

    CREATE TABLE IF NOT EXISTS location
        (
            Id_location INTEGER PRIMARY KEY AUTOINCREMENT,
            region      TEXT,
            region_full TEXT,
            city        TEXT,
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
             Id_contact INTEGER NOT NULL REFERENCES contact(Id_contact)
           );

    -- institution
    CREATE TABLE IF NOT EXISTS  institution
           (
             Id_institution  INTEGER PRIMARY KEY AUTOINCREMENT,
             Name            TEXT    ,
             category        TEXT ,                    -- Institute, Department, University, Laboratory, Center, etc.
             parent          INTEGER REFERENCES institution(Id_institution),
             Id_contact      INTEGER NOT NULL REFERENCES contact(Id_contact)   --  ???
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
             Id_contact    INTEGER NOT NULL REFERENCES contact(Id_contact)  ,
             position      TEXT
           );

    -- genomic_region
    CREATE TABLE IF NOT EXISTS  genomic_region
           (
             Id_gen   INTEGER PRIMARY KEY AUTOINCREMENT,
             Name     TEXT    ,
             Id_taxa  INTEGER NOT NULL REFERENCES taxa(Id_taxa),
             beg      INT,
             end      INT,
             Id_seq   INTEGER   REFERENCES seq(Id_seq)   -- just to map the positions, preferably a bio_seq
           );

    -- source
    CREATE TABLE IF NOT EXISTS   source
           (
             Id_source      INTEGER NOT NULL,     -- original taxa
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
             Id_country_cod3  TEXT REFERENCES countries(iso3) ,
             FAX       TEXT
           );

    /*
    -- country
    CREATE TABLE IF NOT EXISTS  country_cod
           (
             Id_country_cod INTEGER PRIMARY KEY AUTOINCREMENT,
             cod3           TEXT  UNIQUE ,
             cod2           TEXT  UNIQUE ,
             region         TEXT,
             continent      TEXT
           );

    CREATE TABLE IF NOT EXISTS  country
           (
             Id_country      INTEGER PRIMARY KEY AUTOINCREMENT,
             Id_country_cod  INTEGER REFERENCES country_cod(Id_country_cod) ,
             Name            TEXT,
             language        TEXT
           );
    */

    -- taxa
    CREATE TABLE IF NOT EXISTS  taxa
           (
             Id_taxa     INTEGER PRIMARY KEY AUTOINCREMENT,
             Name        TEXT    ,
             vulgar      TEXT ,
             parent      INTEGER REFERENCES taxa(Id_taxa),               -- NULL = Root
             Id_rank     INTEGER NOT NULL REFERENCES taxa_rank(Id_rank),
             NCBI_TaxID  TEXT
           );

    -- taxa_category
    CREATE TABLE IF NOT EXISTS  taxa_rank
           (
             Id_rank       INTEGER PRIMARY KEY AUTOINCREMENT,
             Name          TEXT,
             kingdom       TEXT,
             parent        INTEGER REFERENCES taxa_rank(Id_rank)        -- NULL = Root
             -- NCBI          TEXT                 -- superkingdom, ...
           );



    -- to_excel  VIEW
    CREATE VIEW  files AS SELECT  path, format FROM seq_file;
    CREATE VIEW  to_excel AS SELECT  Seq.Name, Len FROM seq;
    CREATE VIEW  all_frag AS SELECT  Name, Len, pbeg, pend FROM seq, aligned_seq
             ON seq.Id_seq = aligned_seq.Id_part ;
    CREATE VIEW  strains_with_duplicate_isolates AS SELECT strain.Name, count(isolate.Id_strain)
                 FROM strain, isolate  where strain.Id_strain=isolate.Id_strain
                 GROUP BY isolate.Id_strain HAVING count(isolate.Id_strain) >1 ORDER BY -count(isolate.Id_strain)
            /*
             Seq_Name  ,
             parent ,
             NCBI
           )
CREATE VIEW  exl AS SELECT  seq.Name,
                            taxa.Name,
                            seq.Len,
							aligned_seq.pbeg ,
							aligned_seq.pend
					FROM classified_seq, seq, aligned_seq,  taxa
					USING (Id_algseq, aligned_seq.Id_part , Id_taxa)		       ;

*/