    -- http://biosql.org/wiki/Schema_Overview
	-- https://github.com/biosql/biosql/blob/master/sql/biosqldb-sqlite.sql

	-- PRAGMA foreign_keys
	
	CREATE TABLE IF NOT EXISTS  seq_file
                     (
                       Id_file   INTEGER PRIMARY KEY AUTOINCREMENT,
                       path      TEXT    ,    -- UNIQUE   ??
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

/*
    -- seq_file_pos
    -- todo: how to implement? to save space in the DB, and for consistence (keeping only one copy of the sequence, externally )
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
*/

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

    -- isolate_seq   :   all the experimental sequences obtained from a given isolate
    CREATE TABLE IF NOT EXISTS  isolate_seq
           (
             Id_isolate_seq   INTEGER   PRIMARY KEY AUTOINCREMENT,
             Id_isolate       INTEGER   REFERENCES isolate(Id_isolate),
             Id_seq           INTEGER   REFERENCES seq(Id_seq)    -- PRIMARY KEY  ??
           );

    CREATE TABLE IF NOT EXISTS  classified_seq
           (
             Id_clas_seq      INTEGER PRIMARY KEY AUTOINCREMENT,
             Id_taxa          INTEGER   REFERENCES taxa,          -- the finest available classification
             -- description      TEXT,                               --  ??
             -- Id_isolate       INTEGER   REFERENCES isolate    ,
             Id_algseq        INTEGER  NOT NULL REFERENCES aligned_seq(Id_algseq)   -- PRIMARY KEY  ???????
           );

    CREATE TABLE IF NOT EXISTS  ref_schema
           (
             Id_ref_schema    INTEGER PRIMARY KEY AUTOINCREMENT,
             schema           TEXT UNIQUE,        -- PRIMARY KEY  ?   Lu, VRA, IC, etc.
             name             TEXT                -- full text        Lu et. al. 2006, Vina-Rodriguez 2016, etc.
           );

    -- ref_seq  :  name may be use to pre-construct a set of ref seq, and posteriorly updated with actual Id_seq.
    CREATE TABLE IF NOT EXISTS  ref_seq
           (
             Id_ref_seq     INTEGER PRIMARY KEY AUTOINCREMENT,
             Id_taxa        INTEGER   REFERENCES taxa,                    -- the finest available classification
             Id_ref_schema  INTEGER NOT NULL REFERENCES ref_schema(Id_ref_schema),
             Id_seq         INTEGER REFERENCES seq(Id_seq),               -- partial_seq,
             name           TEXT                                          -- preferably = seq.name, but not obligatory
             -- Id_isolate       INTEGER   REFERENCES isolate    ,
             -- Id_algseq        INTEGER  NOT NULL REFERENCES aligned_seq
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
             Id_strain      INTEGER PRIMARY KEY AUTOINCREMENT,
             Name           TEXT,                       --    UNIQUE,      --  ??
             Id_taxa        INTEGER   REFERENCES taxa,       -- the finest available classification, consensus
             host           TEXT,      -- todo: Id_host     INTEGER, --NOT NULL,  -- original taxa    -- consensus
             source         TEXT,      -- todo: Id_source   INTEGER, --NOT NULL,  -- consensus
             year           INT,                             -- consensus
             Id_country_cod TEXT REFERENCES countries(iso3)  -- consensus
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
             country_iso3    TEXT REFERENCES countries(iso3),
             region          TEXT,      -- todo: Id_location,
             region_full     TEXT       -- todo: Id_location
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

    -- non-normalized table, "temporal", just to make possible to search the tree without recursion.
    -- for each inserted taxa here we have all the parents with the rank. Including self !!!
    CREATE TABLE IF NOT EXISTS  taxa_parents
           (
             Id_taxa     INTEGER NOT NULL REFERENCES taxa(Id_taxa),
             parent      INTEGER NOT NULL REFERENCES taxa(Id_taxa),               -- NULL = Root
             Id_rank     INTEGER NOT NULL REFERENCES taxa_rank(Id_rank),
             PRIMARY KEY (Id_taxa, parent)
           );




    -- to_excel  VIEW
    CREATE VIEW  files AS SELECT  path, format FROM seq_file;
    CREATE VIEW  all_frag AS SELECT  Name, Len, pbeg, pend FROM seq, aligned_seq
             ON seq.Id_seq = aligned_seq.Id_part ;
    CREATE VIEW  strains_with_duplicate_isolates AS SELECT strain.Name, count(isolate.Id_strain) AS str_c
                 FROM strain JOIN isolate USING (Id_strain)
                 GROUP BY isolate.Id_strain HAVING count(isolate.Id_strain) >1 ORDER BY str_c DESC;
    CREATE VIEW taxa_tree AS
    SELECT (SELECT Name FROM taxa WHERE Id_taxa=x.Id_taxa ) AS Taxa_name,
           (SELECT Name FROM taxa WHERE Id_taxa=x.parent ) AS Taxa_parent,
           (SELECT Name FROM taxa_rank WHERE Id_rank=x.Id_rank ) AS Taxa_rank
    FROM taxa_parents as x  ORDER BY Id_taxa, Id_rank DESC;
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
    CREATE VIEW  to_excel AS
    SELECT Seq.Name                                         AS MEGA_name,
           t.Name                                           AS Classification,
           (SELECT taxa.Name FROM taxa_parents
                             JOIN taxa_rank  USING (Id_rank)
                             JOIN taxa       ON (taxa.Id_taxa=taxa_parents.parent)
                             WHERE taxa_parents.Id_taxa  =c.Id_taxa
                               AND taxa_rank.Name='genotype')
                                                            AS Genotype,
           (SELECT taxa.Name FROM taxa_parents
                             JOIN taxa_rank  USING (Id_rank)
                             JOIN taxa       ON (taxa.Id_taxa=taxa_parents.parent)
                             WHERE taxa_parents.Id_taxa  =c.Id_taxa
                               AND taxa_rank.Name='group')
                                                            AS 'Group',
           (SELECT taxa.Name FROM taxa_parents
                             JOIN taxa_rank  USING (Id_rank)
                             JOIN taxa       ON (taxa.Id_taxa=taxa_parents.parent)
                             WHERE taxa_parents.Id_taxa  =c.Id_taxa
                               AND taxa_rank.Name='subtype')
                                                            AS Subtype,
           strain.Name                                      AS 'Strain',
           isolate.Name                                     AS 'Isolate ' ,
           countries.name                                   AS 'Country ',
           countries.iso3                                   AS 'Country cod',
           isolate.host                                     AS 'Host',
           isolate.source                                   AS 'Source ',
           isolate.year	                                    AS 'Year',
           isolate.month                                    AS 'Month',
           isolate.day                                      AS 'Day',
           isolate.institution                              AS 'Inst',
           seq.Len                                          AS 'Length ',
           aligned_seq.pbeg                                 AS 'Beg ',
           aligned_seq.pend                                 AS 'End ',
           aligned_seq.pend - aligned_seq.pbeg + 1          AS 'Al length'

    FROM classified_seq AS c JOIN aligned_seq USING (ID_algseq)
                             JOIN seq         ON    (ID_seq=ID_part)
                             JOIN isolate_seq USING (ID_seq)
                             JOIN isolate     USING (ID_isolate)
                             JOIN taxa AS t   USING (Id_taxa)
                             JOIN strain      USING (ID_strain)
                             JOIN countries   ON    (isolate.country_iso3=countries.iso3)
    ;
