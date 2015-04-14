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
	start INTEGER NOT NULL, -- default 0.
	length INTEGER, -- if null, assume sequence.lenght - start 
	md5checksum TEXT, -- if null, assume sequence.md5checksum
	isDerived BOOLEAN NOT NULL, 
	sourceDivergence REAL, -- may be null if isDerived is FALSE (not sure if it needs to be stored)?
	ncbiTaxonID INTEGER, -- ID from http://www.ncbi.nlm.nih.gov/taxonomy, may be null
	isPrimary BOOLEAN NOT NULL,
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

--
-- FIXME: All tables below this point are extremely tentative
--

-- Tables common to both site and allelic models
CREATE TABLE VariantSet (ID INTEGER PRIMARY KEY,
	referenceSetID INTEGER NOT NULL,
	-- datasetID -- TODO: what is this? What table would it point to?
	vcfURI TEXT NOT NULL,
	FOREIGN KEY(referenceSetID) REFERENCES ReferenceSet(ID));
CREATE TABLE CallSet (ID INTEGER PRIMARY KEY,
	name TEXT, -- can be null?
	sampleID TEXT); -- Can the rest of the info now be obtained from VCF?
CREATE TABLE Sequence_VariantSet_Join (sequenceID INTEGER NOT NULL, 
	variantSetID INTEGER NOT NULL,
	PRIMARY KEY(sequenceID,variantSetID),
	FOREIGN KEY(sequenceID) REFERENCES Sequence(ID),
	FOREIGN KEY(variantSetID) REFERENCES VariantSet(ID));
CREATE TABLE VariantSet_CallSet_Join (variantSetID INTEGER NOT NULL, 
	callSetID INTEGER NOT NULL,
	PRIMARY KEY(variantSetID, callSetID),
	FOREIGN KEY(variantSetID) REFERENCES VariantSet(ID),
	FOREIGN KEY(callSetID) REFERENCES CallSet(ID));

-- Site (classical) model tables 
-- We shouldn't be storing anything here besides what's needed for the Allele table join.
CREATE TABLE Variant (ID INTEGER PRIMARY KEY, 
	variantName TEXT NOT NULL, -- corresponding variant in the VCF, accessible via index?
	variantSetID INTEGER NOT NULL,
	FOREIGN KEY(variantSetID) REFERENCES VariantSet(ID)); -- any metadata needed here?
-- is the Call table even needed here? Can't it all be obtained from the VCF?
CREATE TABLE Call (callSetID INTEGER, 
	variantID INTEGER NOT NULL,
	PRIMARY KEY(callSetID, variantID),
	FOREIGN KEY(callSetID) REFERENCES CallSet(ID),
	FOREIGN KEY(variantID) REFERENCES Variant(ID));

-- Allelic model tables
CREATE TABLE Allele (ID INTEGER PRIMARY KEY, 
	variantSetID INTEGER,  -- TODO: Can this be null? cf. Heng Li's question about Alleles and VariantSets
	FOREIGN KEY(variantSetID) REFERENCES VariantSet(ID));
CREATE TABLE AllelePathItem (alleleID INTEGER, 
	pathItemIndex INTEGER NOT NULL, -- one-based index of this pathItem within the entire path
	sequenceID INTEGER NOT NULL, start INTEGER NOT NULL,
	length INTEGER NOT NULL, strandIsForward BOOLEAN NOT NULL,
	PRIMARY KEY(alleleID, pathItemIndex),
	FOREIGN KEY(alleleID) REFERENCES allele(ID),
	FOREIGN KEY(sequenceID) REFERENCES Sequence(ID));
CREATE TABLE Allele_Variant_Join (alleleID INTEGER NOT NULL,
	variantID INTEGER NOT NULL,
	PRIMARY KEY(alleleID, variantID),
	FOREIGN KEY(alleleID) REFERENCES allele(ID),
	FOREIGN KEY(variantID) REFERENCES Variant(ID));
CREATE TABLE AlleleCall (alleleID INTEGER NOT NULL, 
	callSetID INTEGER, variantID INTEGER, -- TODO: Can these be null?
	ploidy INTEGER NOT NULL,
	PRIMARY KEY(alleleID, callSetID, variantID),
	FOREIGN KEY(alleleID) REFERENCES allele(ID),
	FOREIGN KEY(callSetID) REFERENCES CallSet(ID),
	FOREIGN KEY(variantID) REFERENCES Variant(ID)); -- TODO: add metadata? Also, is the PK correct here?


