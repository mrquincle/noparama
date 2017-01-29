
typedef std::vector<data_t> ordered_data_t;

class MCMC {
	private:

		membertrix _membertrix;

		ordered_data_t _dataset;
	public:
		MCMC();

		void run(ordered_data_t & dataset);
};
