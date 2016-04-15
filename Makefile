HOCKING-glm-optimality.pdf: HOCKING-glm-optimality.tex figure-interval-loss-derivative.pdf figure-lasso-criteria.pdf
	pdflatex HOCKING-glm-optimality
figure-interval-loss-derivative.pdf: figure-interval-loss-derivative.R
	R --no-save < $<
figure-lasso-criteria.pdf: figure-lasso-criteria.R
	R --no-save < $<
