/*!
 * Maybe this class doesn't make sense. If we have t_distribution instead with () as function, we might have enough.
 */
class t_sample_pdf {
	private:
		//! hyper parameters
		SufficientStatistics _sufficient_statistics;

	public:

		t_sample_pdf(SufficientStatistics & hyper);

		//
		sample();

};
