// "Cleans" tapes by finding the most likely path through a hidden Markov model (HMM),
// using the Viterbi algorithm. Usually works much worse than e.g. TAPclean;
// you have been warned :-)
//
// Takes in a cycles.plot (from decode) on stdin.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>

#include "tap.h"

#define STATE_BREADTH 75

// From SLC, reportedly from measuring the CIA chip.
#define SYNC_REFERENCE 378.0

// From http://c64tapes.org/dokuwiki/doku.php?id=loaders:rom_loader.
#define ROM_SHORT_PULSE_LENGTH (0x30 * 8)  // 384
#define ROM_MEDIUM_PULSE_LENGTH (0x42 * 8)  // 528
#define ROM_LONG_PULSE_LENGTH (0x56 * 8)  // 688

// From TAPclean sources.
#define NOVA_SHORT_PULSE_LENGTH 288
#define NOVA_LONG_PULSE_LENGTH 688

using namespace std;

enum State { STATE_ROM_SYNC, STATE_ROM_LOADER, STATE_NOVA, STATE_PAUSE };
enum ROMPulseType { ROM_PULSE_TYPE_SHORT = 0, ROM_PULSE_MEDIUM, ROM_PULSE_LONG };
enum NOVAPulseType { NOVA_PULSE_SHORT = 0, NOVA_PULSE_LONG };

static int total_alloc = 0;
		
static const double rom_lengths[] = { ROM_SHORT_PULSE_LENGTH, ROM_MEDIUM_PULSE_LENGTH, ROM_LONG_PULSE_LENGTH };
static const double nova_lengths[] = { NOVA_SHORT_PULSE_LENGTH, NOVA_LONG_PULSE_LENGTH };

struct Node {
	Node() : prev_node(NULL), refcount(1) { ++total_alloc; }
	void ref() const { ++refcount; }
	void unref() const { if (--refcount == 0) { if (prev_node) { prev_node->unref(); } delete this; --total_alloc; } }
	void set_prev_node(const Node *node) { if (prev_node) { prev_node->unref(); } node->ref(); prev_node = node; }
	const Node *get_prev_node() const { return prev_node; }

	// Viterbi information.
	double cost;
private:
	const Node* prev_node;
public:
	int emit;
	
	// State information.
	State state;

	union {
		// For STATE_ROM_SYNC only.
		struct {
			int num_pulses_left;
		};

		// For STATE_ROM_LOADER.
		struct {
			ROMPulseType last_pulse_type;
			bool first_in_pair;  // If the _last_ pulse was first in the pair.
		};

		// For STATE_NOVA.
		struct {
			NOVAPulseType last_nova_pulse_type;
		};
	};

private:
	mutable int refcount;
	~Node() {}
	Node(const Node &);
};

struct StateLessThanComparator {
	bool operator() (const Node *a, const Node *b) const
	{
		if (a->state != b->state)
			return (a->state < b->state);

		if (a->state == STATE_ROM_SYNC) {
			if (a->num_pulses_left != b->num_pulses_left)
				return (a->num_pulses_left < b->num_pulses_left);
		}

		if (a->state == STATE_ROM_LOADER) {
			if (a->last_pulse_type != b->last_pulse_type)
				return (a->last_pulse_type < b->last_pulse_type);
			if (a->first_in_pair != b->first_in_pair)
				return (a->first_in_pair < b->first_in_pair);
		}
		
		if (a->state == STATE_NOVA) {
			if (a->last_nova_pulse_type != b->last_nova_pulse_type)
				return (a->last_nova_pulse_type < b->last_nova_pulse_type);
		}

		return false;
	}
};

struct StateEqualsComparator {
	StateLessThanComparator lt;

	bool operator() (const Node *a, const Node *b) const
	{
		return !(lt(a, b) || lt(b, a));
	}
};

struct CompareByCost {
	bool operator() (const Node *a, const Node *b) const
	{
		return a->cost < b->cost;
	}
};

double penalty(double length, double reference, double ratio, double reference_ratio)
{
	double ratio_penalty = ((ratio > reference_ratio) ? ratio / reference_ratio : reference_ratio / ratio) - 1.0;
	double distance_penalty = fabs(length - reference);
	return ratio_penalty * ratio_penalty + 1e-4 * distance_penalty;
}

double pulse_penalty(int num_pulses_left)
{
	int num_pulses = 27136 - num_pulses_left;
	return fabs(sqrt(27136.0) - sqrt(num_pulses));
	// min(abs(prev->num_pulses_left - (27136 - 0x4f)), abs(prev->num_pulses_left));
}

void possibly_add_state(Node *node, vector<Node *> *next_states)
{
	assert(node->cost >= node->get_prev_node()->cost);
	next_states->push_back(node);
}
		
void extend_state(const Node *prev, double last_length, double length, double ratio, vector<Node *> *next_states)
{
	// If this pulse is  long, it means we could transition into another type.
	if (length > 700.0) {
		double cost = 0.0;
		if (prev->state == STATE_ROM_SYNC) {
			// We don't really want to jump directly from sync to Novaload sync...
			cost = pulse_penalty(prev->num_pulses_left);
		} else if (prev->state == STATE_ROM_LOADER) {
			// Jumping in the middle of a bit is bad, too
			cost = 50.0;
		} else if (prev->state == STATE_NOVA) {
			// is this too close to a long bit?
			double prev_ref = nova_lengths[prev->last_nova_pulse_type];
			double long_pulse_would_be = last_length * (NOVA_LONG_PULSE_LENGTH / prev_ref);
			double ratio_penalty = max(2.5 - length / long_pulse_would_be, 0.0);
			double distance_penalty = max(1000.0 - length, 0.0);
			cost = ratio_penalty * ratio_penalty + 1e-4 * distance_penalty;
		} else {
			// Make sure that marking things as pauses are not _that_ painless...
			cost = 1e-4 * max(1000.0 - length, 0.0);
		}

		Node *n = new Node;
		n->cost = prev->cost + cost; // + 10.0;
		n->set_prev_node(prev);
		n->emit = length;

		n->state = STATE_PAUSE;
		possibly_add_state(n, next_states);
	}

	// If in STATE_ROM_SYNC, it's possible that this was another sync pulse.
	if (prev->state == STATE_ROM_SYNC || prev->state == STATE_PAUSE) {
		Node *n = new Node;
		if (prev->state == STATE_PAUSE) {
			n->cost = prev->cost + penalty(length, SYNC_REFERENCE, 1.0, 1.0);
		} else {
			n->cost = prev->cost + penalty(length, SYNC_REFERENCE, ratio, 1.0);
		}
		n->set_prev_node(prev);
		n->emit = SYNC_REFERENCE;

		n->state = STATE_ROM_SYNC;	
		if (prev->state == STATE_PAUSE) {
			n->num_pulses_left = 27136;
		} else {
			n->num_pulses_left = prev->num_pulses_left - 1;
		}
		possibly_add_state(n, next_states);
	}

	// If in STATE_ROM_SYNC, maybe we transitioned into ROM_LOADER.
	// That always starts with a (L,M).
	if (prev->state == STATE_ROM_SYNC) {
		Node *n = new Node;
		n->cost = prev->cost +
			penalty(length, ROM_LONG_PULSE_LENGTH, ratio, ROM_LONG_PULSE_LENGTH / SYNC_REFERENCE) +
			0.1 * pulse_penalty(prev->num_pulses_left);
		n->set_prev_node(prev);
		n->emit = ROM_LONG_PULSE_LENGTH;

		n->state = STATE_ROM_LOADER;
		n->last_pulse_type = ROM_PULSE_LONG;
		n->first_in_pair = true;
		possibly_add_state(n, next_states);
	}
	
	// If in ROM_LOADER, we could have short, medium or long pulses.
	if (prev->state == STATE_ROM_LOADER) {
		for (int pulse_type = ROM_PULSE_TYPE_SHORT; pulse_type <= ROM_PULSE_LONG; ++pulse_type) {

			// Filter illegal ROM loader pairs.
			if (prev->last_pulse_type == ROM_PULSE_LONG && pulse_type == ROM_PULSE_LONG) {
				continue;
			}
			if (prev->first_in_pair) {
				if (prev->last_pulse_type == ROM_PULSE_TYPE_SHORT && pulse_type != ROM_PULSE_MEDIUM) {
					continue;
				}
				if (prev->last_pulse_type == ROM_PULSE_MEDIUM && pulse_type != ROM_PULSE_TYPE_SHORT) {
					continue;
				}
				if (pulse_type == ROM_PULSE_LONG) {
					continue;
				}
			}

			Node *n = new Node;
			n->cost = prev->cost +
				penalty(length, rom_lengths[pulse_type], ratio, rom_lengths[pulse_type] / rom_lengths[prev->last_pulse_type]);
			n->set_prev_node(prev);
			n->emit = rom_lengths[pulse_type];

			if (prev->first_in_pair && (prev->last_pulse_type == ROM_PULSE_LONG && pulse_type == ROM_PULSE_TYPE_SHORT)) {
				// (L,S) = end-of-data-marker
				n->state = STATE_ROM_SYNC;
				n->num_pulses_left = 0x4f;  // http://c64tapes.org/dokuwiki/doku.php?id=loaders:rom_loader
			} else {
				n->state = STATE_ROM_LOADER;
				n->last_pulse_type = ROMPulseType(pulse_type);
				n->first_in_pair = !prev->first_in_pair;
			}
			possibly_add_state(n, next_states);
		}		
	}

	// If in STATE_NOVA, we only have long and short pulses.
	if (prev->state == STATE_NOVA || prev->state == STATE_PAUSE || prev->state == STATE_ROM_SYNC) {
		for (int pulse_type = NOVA_PULSE_SHORT; pulse_type <= NOVA_PULSE_LONG; ++pulse_type) {
			Node *n = new Node;
			if (prev->state == STATE_PAUSE) {
				// going from PAUSE to NOVA is only realistic after a pretty long pause.
				n->cost = prev->cost + fabs(2000.0f - last_length) + penalty(length, nova_lengths[pulse_type], 1.0, 1.0);
			} else if (prev->state == STATE_ROM_SYNC) {
				n->cost = prev->cost +
					penalty(length, NOVA_SHORT_PULSE_LENGTH, ratio, nova_lengths[pulse_type] / SYNC_REFERENCE) +
					0.1 * pulse_penalty(prev->num_pulses_left);
			} else {
				n->cost = prev->cost +
					penalty(length, nova_lengths[pulse_type], ratio, nova_lengths[pulse_type] / nova_lengths[prev->last_nova_pulse_type]);
			}
			n->set_prev_node(prev);
			n->emit = nova_lengths[pulse_type];

			n->state = STATE_NOVA;
			n->last_nova_pulse_type = NOVAPulseType(pulse_type);
			possibly_add_state(n, next_states);
		}
	}

	// TODO: Other loader types
}

int main(int argc, char **argv)
{
	Node *start = new Node;
	start->cost = 0.0;
	start->emit = -1;
	start->state = STATE_PAUSE;

	vector<Node *> states;
	states.push_back(start);

	double last_length = SYNC_REFERENCE;

	int pulse_num = 0;
	int max_total_alloc = 0;
	vector<pair<double, double>> original_pulses;
	for ( ;; ) {
		double time, length;
		if (scanf("%lf %lf", &time, &length) != 2) {
			break;
		}
		original_pulses.push_back(make_pair(time, length));

		++pulse_num;

		max_total_alloc = max(total_alloc, max_total_alloc);
		if (total_alloc > 20000000) {
			printf("More than 20M states reached (out of RAM); aborting at pulse %d.\n", pulse_num);
			break;
		}

		//if (length > 2000) {
		//	break;
		//}

		double ratio = length / last_length;
	
		vector<Node *> next_states;
		for (unsigned i = 0; i < states.size(); ++i) {
			extend_state(states[i], last_length, length, ratio, &next_states);

			// We no longer need this state; if it doesn't have any children,
			// it can go away.
			states[i]->unref();
		}
		states.clear();

		// Remove duplicates, tie-breaking by score.
		sort(next_states.begin(), next_states.end(), CompareByCost());
		stable_sort(next_states.begin(), next_states.end(), StateLessThanComparator());

		// unique and move in one step. Do not use std::unique(),
		// it has very wrong behavior for pointers!
		for (unsigned i = 0; i < next_states.size(); ++i) {
			if (i > 0 && StateEqualsComparator()(next_states[i], states.back())) {
				next_states[i]->unref();
			} else {
				states.push_back(next_states[i]);
			}
		}
		assert(!states.empty());
			
		// Prune unlikely next_states to save time and memory.
		if (states.size() >= STATE_BREADTH) {
			sort(states.begin(), states.end(), CompareByCost());  // nth element?
			for (unsigned i = STATE_BREADTH; i < states.size(); ++i) {
				states[i]->unref();
			}
			states.resize(STATE_BREADTH);
		}
		last_length = length;

		if (pulse_num % 1000 == 0) {
			fprintf(stderr, "\nProcessing pulses... %d [%5.2f,%7.2f]  ", pulse_num, time, length);
			for (unsigned i = 0; i < 3 && i < states.size(); ++i) {
				fprintf(stderr, " [%d:state=%d cost=%f]", i, states[i]->state, states[i]->cost);
			}
		}	
	}

	// Find the best final node.
	const Node *best_node = NULL;
	for (unsigned i = 0; i < states.size(); ++i) {
		if (states[i]->state == STATE_ROM_SYNC) {
			states[i]->cost += 0.1 * pulse_penalty(states[i]->num_pulses_left);
		}

		if (best_node == NULL || states[i]->cost < best_node->cost) {
			best_node = states[i];
		}
	}

	fprintf(stderr, "\rTotal cost is %f. Peak RAM usage %.2f MB, pluss malloc overhead.\n",
		best_node->cost, sizeof(Node) * max_total_alloc / 1048576.0);

	// Backtrack.
	vector<const Node *> cleaned;
	while (best_node != NULL) {
		cleaned.push_back(best_node);
		best_node = best_node->get_prev_node();
	}

	reverse(cleaned.begin(), cleaned.end());

	for (unsigned i = 0; i < cleaned.size(); ++i) {
		//fprintf(stderr, "%d: state %d emit %d cost %f\n", i, cleaned[i]->state, cleaned[i]->emit, cleaned[i]->cost);
		if (cleaned[i]->emit <= 0) {
			continue;
		}
		int state;
		if (cleaned[i]->state == STATE_ROM_SYNC) {
			state = 0;
		} else if (cleaned[i]->state == STATE_ROM_LOADER) {
			state = 1 + cleaned[i]->last_pulse_type;
		} else if (cleaned[i]->state == STATE_NOVA) {
			state = 4 + cleaned[i]->last_nova_pulse_type;
		} else if (cleaned[i]->state == STATE_PAUSE) {
			state = 6;
		} else {
			assert(false);
		}
	//	printf("%d\n", cleaned[i]->emit);
		printf("%f %f %d\n", original_pulses[i - 1].first, original_pulses[i - 1].second, state);
	}

	// output TAP file	
	FILE *fp = fopen("cleaned.tap", "wb");
	std::vector<char> tap_data;
	for (unsigned i = 0; i < cleaned.size(); ++i) {
		if (cleaned[i]->emit <= 0) {
			continue;
		}
		int len = lrintf(cleaned[i]->emit / TAP_RESOLUTION);
		if (len <= 255) {
			tap_data.push_back(len);
		} else {
			int overflow_len = lrintf(cleaned[i]->emit);
			tap_data.push_back(0);
			tap_data.push_back(overflow_len & 0xff);
			tap_data.push_back((overflow_len >> 8) & 0xff);
			tap_data.push_back(overflow_len >> 16);
		}
	}

	tap_header hdr;
	memcpy(hdr.identifier, "C64-TAPE-RAW", 12);
	hdr.version = 1;
	hdr.reserved[0] = hdr.reserved[1] = hdr.reserved[2] = 0;
	hdr.data_len = tap_data.size();

	fwrite(&hdr, sizeof(hdr), 1, fp);
	fwrite(tap_data.data(), tap_data.size(), 1, fp);
	fclose(fp);
}
