--
-- Add extra terms used by Artemis with Chado
--

-- 
-- add term for evidence codes
--
INSERT INTO dbxref 
     ( db_id, accession ) 
VALUES
     ( (SELECT db_id FROM db WHERE name = 'null'), 'Typically an evidence code' );

INSERT INTO cvterm 
     ( cv_id,  name, dbxref_id ) 
VALUES 
     ( (SELECT cv_id FROM cv WHERE name ='feature_property'), 'evidence', 
       (SELECT dbxref_id FROM dbxref WHERE accession='Typically an evidence code') );

--
-- database used for controlled curation terms
--
INSERT INTO db ( name ) VALUES (  'CCGEN' );

-- 
-- add terms for literature
--
INSERT INTO cv
 ( name, definition) 
VALUES
 ( 'genedb_literature', 'terms for literature' );

INSERT INTO dbxref 
 ( db_id, accession ) 
VALUES  
 ( (SELECT db_id FROM db WHERE name = 'null'), 'unfetched literature type');

INSERT INTO dbxref 
 ( db_id, accession ) 
VALUES  
 ( (SELECT db_id FROM db WHERE name = 'null'), 'journal literature type' );

INSERT INTO dbxref 
 ( db_id, accession ) 
VALUES   
 ( (SELECT db_id FROM db WHERE name = 'null'), 'unknown literature type' );

INSERT INTO cvterm 
 ( cv_id, name, dbxref_id ) 
VALUES 
 ( (SELECT cv_id FROM cv WHERE name = 'genedb_literature'), 'unfetched', 
   (SELECT dbxref_id FROM dbxref WHERE accession='unfetched literature type') );

INSERT INTO cvterm 
 ( cv_id, name, dbxref_id ) 
VALUES 
 ( (SELECT cv_id FROM cv WHERE name = 'genedb_literature'), 'journal', 
   (SELECT dbxref_id FROM dbxref WHERE accession='journal literature type') );
       
INSERT INTO cvterm 
 ( cv_id, name, dbxref_id ) 
VALUES 
 ( (SELECT cv_id FROM cv WHERE name = 'genedb_literature'), 'unknown', 
   (SELECT dbxref_id FROM dbxref WHERE accession='unknown literature type') );


