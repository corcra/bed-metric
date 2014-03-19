bed-metric
==========

Calculate a sort of 'information content' for data in bed format, by seeing how many unique positions are recovered as a function of the total read depth.

==========
Usage: python bed-metric.py FILE.gz (--plot)

Optional argument --plot automatically generates a pdf of the unique positions
as a function of total positions, with 'predicted asymptotic' value if
necessary.
