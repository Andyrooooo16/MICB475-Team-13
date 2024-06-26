# Meeting Notes 2024/03/14

- Order x-labels (from low to high)
- Consistent across metrics in terms of alpha diversity
- Potential explain trends even if p-value isn’t significant (ideally if p-value is close to significant)
	○ Clearly state non-significant data if going onto expalin trends
- Pair wise comparisons b/w two key groups 
	○ Inflammation groups w/ and w/o surgery
	○ Surgery w/ and w/o inflammation
- Calprotectin vs. Disease severity
	○ Switch calprotectin to inflammation
- Create a new metadata column with disease severity and presence or absence of inflammation
	○ Can either do in excel
	○ Or before phyloseq object in R
- Red for inflammation 
- Blue for no inflammation 
- Strongest colour for high disease severity
- Relative abundance:
	○ Group smaller phyla together (<1%)
		§ Follow up email on how to do this
- Taxa bar plot:
	○ Note observations for paper (correlation not causation)
- Indicator species: 
	○ Avril will try and update us 
- Core microbiome
	○ Only include taxa from indicator species (key species) and find where they are using cd_locaiton
	○ Make sure to explain in manuscript how this aim connects to other aims.

- Remove failed samples from phyloseq
	○ Mention why they were removed in manuscript
