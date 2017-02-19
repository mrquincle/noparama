#include <string>
#include <membertrix>

typedef std::vector<int> ground_truth_t;

class Results {
    private:
        const membertrix & _membertrix;

        ground_truth_t & _ground_truth;

        int _verbosity;
    public:
        Results(const membertrix &membertrix, ground_truth_t & ground_truth);

        void write(const std::string & path, const std::string & basename);
};
