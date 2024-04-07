# Meeting notes 2024/03/28

**Meeting 6 Figures for Avril Slide Deck**: 
*Slide 2: Colonic CD*
- Add stars for p-value
- Add Healthy Control

*Slide 3: Groups with low disease severity have higher alpha diversity indexes*
- Add Health control to boxplot if possible

*Slide 6: Inflammation has no effect on alpha diversity*
- Facet
- Add p-values only for significant groups
- bubble plot for ISA --> requires some data wrangling (need to calculate relative abundance ourselves, take average abundance for each group and species) (should be able to do in 2 lines of code) (use geom point) (use groupby and summarize functions
.7 and above for the cutoff

*Slide 7: Bubble Plot*
- Fix labels, add axes, add legend
- ISA performed on inflammation w/ surgery
- Used relative abundance to generate bubble plot
- First three indicator species, add family name instead of ‘uncultured’ → need to manually add names
- Is this the best way to show data? 
- Bubble plots should be coloured by category (if indicator species associated with category, assign it a colour; bubble size = RA)
- Typically shows species that is only in one group
- Initially ISA analysis returned 48, when using RA of 0.8 or higher, returned 18 
- Taxa that stand out and seem to belong to just one category:
- Christensenellaceae
- Agathobacter (maybe)
- UCG-002 → Need to identify family

*Slide 9: Core Microbiome - Inflammation and Surgery*
- How to change colour scheme? How to convert % to species/genera? 



