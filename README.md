# ComputerScience478470
Computer Science Course Paper

This project is for the Course Computer Science for Business Analytics of the Erasmus University Rotterdam. The challenge was to find duplicates in a set of Tv's where scalability was an important part. Therefor, Locality Sensitivity Hashing needs to be used. 

The code starts with data cleaning, after which the Model Word set is constructed. After this part, the binary matrix is filled in for all Tv's, giving a 1 if a model word is present in the title of that Tv.
Then the binary matrix is hashed to a signature matrix which is used as input in the LSH algorithm which presents candidate pairs.
Those candidate pairs are then compared using a MSM algorithm. The F1* and F1 are reported. 

The code can be used on the dataset "TVs-all-merged.json" without entering anything else. 
