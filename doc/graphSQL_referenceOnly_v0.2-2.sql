-- Sequence table and FASTA URI reference - basis of graph topology
CREATE TABLE FASTA (ID INTEGER PRIMARY KEY,
	fastaURI TEXT NOT NULL);
CREATE TABLE Sequence (ID INTEGER PRIMARY KEY,
	-- FASTA file URI and record name within it containing sequence  
	fastaID INTEGER, -- If null, query enclosing ReferenceSet's fastaID
	sequenceRecordName TEXT NOT NULL, 
	md5checksum TEXT NOT NULL,
	length INTEGER NOT NULL,
	FOREIGN KEY(fastaID) REFERENCES FASTA(ID));
CREATE TABLE GraphJoin (ID INTEGER PRIMARY KEY,
	-- startJoin defines parent sequence & position that 5' end attaches to.
	side1SequenceID INTEGER NOT NULL,
	side1Position INTEGER NOT NULL,
	side1StrandIsForward BOOLEAN NOT NULL,
	-- endJoin defines parent sequence & position that 3' end attaches to.
	side2SequenceID INTEGER NOT NULL,
	side2Position INTEGER NOT NULL,
	side2StrandIsForward BOOLEAN NOT NULL,
	FOREIGN KEY(side1SequenceID) REFERENCES Sequence(ID),
	FOREIGN KEY(side2SequenceID) REFERENCES Sequence(ID)); 

-- References
CREATE TABLE Reference (ID INTEGER PRIMARY KEY, 
	name TEXT NOT NULL, 
	updateTime DATE NOT NULL, 
	sequenceID INTEGER NOT NULL,
	isDerived BOOLEAN NOT NULL, 
	sourceDivergence REAL, -- may be null if isDerived is FALSE (not sure if it needs to be stored)?
	ncbiTaxonID INTEGER, -- ID from http://www.ncbi.nlm.nih.gov/taxonomy, may be null
	FOREIGN KEY(sequenceID) REFERENCES Sequence(ID)); 
CREATE TABLE ReferenceAccession (ID INTEGER PRIMARY KEY,
	referenceID INTEGER NOT NULL,
	accessionID TEXT NOT NULL,
	FOREIGN KEY(referenceID) REFERENCES Reference(ID)); 

-- Reference sets
CREATE TABLE ReferenceSet (ID INTEGER PRIMARY KEY,
	ncbiTaxonID INT, -- may be null, may differ from ncbiTaxonID of contained Reference record
	description TEXT, -- may be null, but shouldn't be?
	fastaID INTEGER, -- may be null. TODO: What if both this and a member Reference's Sequence fastaID are null?
	assemblyID TEXT, -- may be null, but REALLY shouldn't be?
	isDerived BOOLEAN NOT NULL,
	FOREIGN KEY(fastaID) REFERENCES FASTA(ID));
CREATE TABLE ReferenceSetAccession (ID INTEGER PRIMARY KEY,
	referenceSetID INTEGER NOT NULL,
	accessionID TEXT NOT NULL,
	FOREIGN KEY(referenceSetID) REFERENCES ReferenceSet(ID)); 

CREATE TABLE Reference_ReferenceSet_Join (referenceID INTEGER NOT NULL, 
	referenceSetID INTEGER NOT NULL,
	PRIMARY KEY(referenceID, referenceSetID),
	FOREIGN KEY(referenceID) REFERENCES Reference(ID),
	FOREIGN KEY(referenceSetID) REFERENCES ReferenceSet(ID));

-- Example graph
-- Imagine a FASTA file, full of data!
INSERT INTO FASTA VALUES (1, 'URI://FASTA_FILE_TODO');
-- Define some sequences
INSERT INTO Sequence VALUES (1, 1, 'ref1', 'md5', 240);
INSERT INTO Sequence VALUES (2, 1, 'alt1', 'md5', 34);
-- And some graph topology
INSERT INTO GraphJoin VALUES (1, 1, 50, 'FALSE', 2, 1, 'FALSE');
INSERT INTO GraphJoin VALUES (2, 2, 34, 'FALSE', 1, 60, 'FALSE');

-- Make a ReferenceSet. 
INSERT INTO ReferenceSet VALUES (1, NULL, NULL, NULL, NULL, 'FALSE');
-- Make a couple References
INSERT INTO Reference VALUES (1, 'base reference', date('now'), 1, 'FALSE', NULL, NULL);
INSERT INTO Reference VALUES (2, 'alternate allele', date('now'), 2, 'FALSE', NULL, NULL);
-- Add them to the set
INSERT INTO Reference_ReferenceSet_Join VALUES (1, 1);
INSERT INTO Reference_ReferenceSet_Join VALUES (2, 1);


