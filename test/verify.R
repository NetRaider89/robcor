X = read.table('test/R_correlations.txt')
Y = read.table('test/python_correlations.txt')

n = length(X)

# Compute mean sum of squared errors
msse = norm(X-Y, type='2')^2 / n

cat('Mean Sum of Squared Errors: ', msse, sep='', '\n')
