#include <np_sufficient_statistics.h>

class PosteriorPredictive {
	private:

		SufficientStatistics _sufficient_statistics;

	public:
		PosteriorPredictive(SufficientStatistics &sufficient_statistics);

		~PosteriorPredictive() {
		}

		SufficientStatistics & Get() {
			return _sufficient_statistics;
		}

		// Get probability
		double Probability();
}

