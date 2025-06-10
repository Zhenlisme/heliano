## 2024-04-04
# Main modefication:
Allow users to provide a curated LTS-RTS file to help identify nonautonomous HLEs in another species whose autonomous counterparts do not exist.
## 2023-10-04
Add the '-flank_sim' parameter which allow users to set the cut-off to define false positive LTS/RTS. The lower the value, the more strigent. This value was set to 0.7 in previous versions but it is now set as 0.5 by default.
## 2025-04-17
Resove the The hmmsearch error;
Add that "--table" parameter that allows users to adjust the genetic code of test organisms;
Optimize the LTS/RTS selection. When there are alternative terminal sequences on the same autonomous locus, try to use the one with higher blastn score.
