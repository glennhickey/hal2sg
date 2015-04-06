-- We have our FASTA file right next to this file, so we can reference it by
-- relative path. Probably.
INSERT INTO FASTA VALUES (1, 'file://sequence.fa');

-- We define all the sequences in that file. Get hashes by removing all line
-- breaks and whitespace at <http://removelinebreaks.net/>, then taking the md5
-- with <http://onlinemd5.com/>, and finally converting the hash to lower case
-- with <http://convertcase.net/>
INSERT INTO Sequence VALUES (1, 1, 'chr1', 'f28788f2203d71311b5d8fe81fe1bd2e', 1000);
INSERT INTO Sequence VALUES (2, 1, 'chr1snp1', '7fc56270e7a70fa81a5935b72eacbe29', 1);
INSERT INTO Sequence VALUES (3, 1, 'chr2', '2895d94c7491966ef9df7af5ecf77e9f', 1000);
-- IDs won't always be contiguous and in order
INSERT INTO Sequence VALUES (8, 1, 'chr2ins1', '4833b4fa1627b1ee25f83698f768f997', 30);
INSERT INTO Sequence VALUES (5, 1, 'chr3', '52ae3ef016c60c5e978306e8d3334cd8', 1000);

-- And some graph topology
-- Attach the SNP on chr1 at base 80
INSERT INTO GraphJoin VALUES (1, 1, 79, 'FALSE', 2, 0, 'TRUE');
INSERT INTO GraphJoin VALUES (2, 1, 81, 'TRUE', 2, 0, 'FALSE');
-- Add a tandem duplication around it.
INSERT INTO GraphJoin VALUES (3, 1, 59, 'TRUE', 1, 74, 'FALSE');
-- Add a deletion on chr2
INSERT INTO GraphJoin VALUES (4, 3, 33, 'FALSE', 3, 108, 'TRUE');
-- Add a translocation from chr1 to chr2, which is just a pure insertion in chr2
INSERT INTO GraphJoin VALUES (5, 1, 233, 'TRUE', 3, 310, 'FALSE');
INSERT INTO GraphJoin VALUES (6, 1, 289, 'FALSE', 3, 311, 'TRUE');
-- Add a 30bp point insert on chr2 at position 300
INSERT INTO GraphJoin VALUES (7, 3, 300, 'FALSE', 8, 0, 'TRUE');
INSERT INTO GraphJoin VALUES (8, 3, 301, 'TRUE', 8, 29, 'FALSE');
-- Nothing needs to touch chr3.

-- Make a ReferenceSet. 
INSERT INTO ReferenceSet VALUES (1, NULL, NULL, 1, 'normal', 'FALSE');
-- Make references, but wth IDs in a different order than the sequences.
INSERT INTO Reference VALUES (1, 'chr1', date('now'), 1, 'FALSE', NULL, NULL);
INSERT INTO Reference VALUES (2, 'chr2', date('now'), 3, 'FALSE', NULL, NULL);
INSERT INTO Reference VALUES (3, 'chr3', date('now'), 5, 'FALSE', NULL, NULL);
INSERT INTO Reference VALUES (4, 'chr1snp1', date('now'), 2, 'FALSE', NULL, NULL);
INSERT INTO Reference VALUES (5, 'chr2ins1', date('now'), 8, 'FALSE', NULL, NULL);
-- Add them to the set
INSERT INTO Reference_ReferenceSet_Join VALUES (1, 1);
INSERT INTO Reference_ReferenceSet_Join VALUES (2, 1);
INSERT INTO Reference_ReferenceSet_Join VALUES (3, 1);
INSERT INTO Reference_ReferenceSet_Join VALUES (4, 1);
INSERT INTO Reference_ReferenceSet_Join VALUES (5, 1);

