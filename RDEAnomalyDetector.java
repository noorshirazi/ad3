import java.math.BigInteger;
import java.util.Arrays;

/**
 * Detects anomalies in a series of vectors, using a Recursive Density
 * Estimation technique.
 */
public class RDEAnomalyDetector {
    /**
     * Specifies the number of dimensions that all submitted vectors
     * must have.
     */
    private final int dimensions;

    /**
     * Specifies the number of consecutively submitted vectors that
     * must keep the density below the mean density while no anomaly
     * is happening before an anomaly is deemed to have started.
     */
    private final int anomalyStartLag;

    /**
     * Specifies the number of consecutively submitted vectors that
     * must keep the density at or above the mean density during an
     * anomaly before the anomaly is deemed to have stopped.
     */
    private final int anomalyEndLag;

    /**
     * Invoked when an anomaly starts.
     */
    private final Runnable onAnomalyStart;

    /**
     * Invoked when an anomaly ends.
     */
    private final Runnable onAnomalyEnd;

    /**
     * Records the mean of all submitted vectors.  This field is
     * {@code null} initially, and is used to detect when the first
     * vector is submitted, and handled differently to later vectors.
     */
    private double[] mean;

    /**
     * Records the scalar product of all submitted vectors.
     */
    private double product;

    /**
     * Records the density of all submitted vectors.
     */
    private double density;

    /**
     * Records the mean of all densities computed at each submission
     * of a vector.
     */
    private double meanDensity;

    /**
     * Records how many vectors have been submitted.
     */
    private long contribution;

    /**
     * Records how many vectors have been submitted since the last
     * change of anomalous state.
     */
    private long stability;

    /**
     * Records whether the analysis is reporting an anomalous state.
     */
    private boolean anomaly;

    /**
     * Records whether the last submitted vector placed the density
     * under the mean density.
     */
    private boolean underMean;

    /**
     * Records how many vector submissions have passed since the
     * density crossed over the mean density.
     */
    private long meanSideCounter;

    /**
     * Create an anomaly detector.  The detector will initially be in
     * a normal (non-anomalous) state.
     *
     * @param dimensions the number of dimensions in each submitted vector
     *
     * @param anomalyStartLag the number of consecutively submitted
     * vectors that must keep the density below the mean density while
     * no anomaly is happening before an anomaly is deemed to have
     * started
     *
     * @param anomalyEndLag the number of consecutively submitted
     * vectors that must keep the density at or above the mean density
     * during an anomaly before the anomaly is deemed to have stopped
     *
     * @param onAnomalyStart invoked when an anomaly starts
     *
     * @param onAnomalyEnd invoked when an anomaly ends
     */
    public RDEAnomalyDetector(int dimensions,
			      int anomalyStartLag,
			      int anomalyEndLag,
			      Runnable onAnomalyStart,
			      Runnable onAnomalyEnd) {
	this.dimensions = dimensions;
	this.anomalyStartLag = anomalyStartLag;
	this.anomalyEndLag = anomalyEndLag;
	this.onAnomalyStart = onAnomalyStart;
	this.onAnomalyEnd = onAnomalyEnd;
    }

    private static double sumOfSquares(double[] point) {
	double result = 0.0;
	for (double var : point)
	    result += var * var;
	return result;
    }

    private double sumOfSquares(double[] point, double[] base) {
	double result = 0.0;
	for (int i = 0; i < dimensions; i++) {
	    final double var = point[i] - base[i];
	    result += var * var;
	}
	return result;
    }

    /**
     * Submit a vector.  It must have the same number of dimensions as
     * specified by the constructor.
     *
     * @param point the submitted vector
     *
     * @throws IllegalArgumentException if the vector has the wrong
     * number of parameters
     */
    public void submit(double[] point) {
	if (point.length != dimensions)
	    throw new IllegalArgumentException("vector size: expected " +
					       dimensions + "; received " +
					       point.length);
	if (mean == null) {
	    mean = Arrays.copyOf(point, dimensions);
	    product = sumOfSquares(point);
	    contribution = 1L;
	    stability = 1L;
	    density = 1.0;
	    meanDensity = density;
	    underMean = false;
	    meanSideCounter = 1L;
	} else {
	    /* Update counters. */
	    final long lastContribution = contribution;
	    contribution++;
	    final long lastStability = stability;
	    stability++;

	    /* Update the product. */
	    final double lastProduct = product;
	    product = (lastProduct * lastContribution + sumOfSquares(point))
		/ contribution;

	    /* Update the mean. */
	    final double[] lastMean = Arrays.copyOf(mean, dimensions);
	    for (int i = 0; i < dimensions; i++)
		mean[i] = (lastMean[i] * lastContribution + point[i])
		    / contribution;

	    /* Update density. */
	    final double lastDensity = density;
	    density = 1.0 / (1.0 + sumOfSquares(point, mean)
			     + product - sumOfSquares(mean));

	    /* Compute change in density. */
	    final double deltaDensity = Math.abs(density - lastDensity);

	    /* Compute mean of density. */
	    final double lastMeanDensity = meanDensity;
	    meanDensity = lastStability * lastMeanDensity + density;
	    meanDensity /= stability;
	    meanDensity *= 1.0 - deltaDensity;
	    meanDensity += deltaDensity * density;

	    /* Update the flag recording whether the density is under
	     * the mean. */
	    final boolean lastUnderMean = underMean;
	    underMean = density < meanDensity;

	    /* Check whether we've crossed the threshold, or record
	     * how long it is since we last crossed it. */
	    if (underMean != lastUnderMean) {
		meanSideCounter = 0;
	    } else {
		meanSideCounter++;
	    }

	    /* Now decide whether we need to change state. */
	    if (!anomaly && underMean &&
		meanSideCounter >= anomalyStartLag) {
		anomaly = true;
		stability = 0;
		onAnomalyStart.run();
	    } else if (anomaly && !underMean &&
		       meanSideCounter >= anomalyEndLag) {
		anomaly = false;
		stability = 0;
		onAnomalyEnd.run();
	    }
	}
    }

    /**
     * Determine whether the detector is in an anomalous state.
     *
     * @return {@code true} if the detector is in an anomalous state;
     * {@code false} otherwise
     */
    public boolean isAnomalous() { return anomaly; }
}
