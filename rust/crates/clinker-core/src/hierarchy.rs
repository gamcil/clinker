//! Agglomerative hierarchy construction and leaf ordering for cluster ordering.

/// Build a Ward hierarchy, then approximately minimize neighbouring distances
/// while preserving its groups. Retaining 32 endpoint alternatives per subtree
/// bounds the leaf-ordering work; the hierarchy uses O(n²) time and memory.
pub(crate) fn ward_leaf_order(distances: &[f64], count: usize) -> Vec<usize> {
    if count < 2 {
        return (0..count).collect();
    }
    bounded_leaf_order(&ward_tree(distances, count), distances, count, 32)
}

fn ward_tree(distances: &[f64], count: usize) -> Node {
    assert_eq!(distances.len(), count * (count - 1) / 2);

    // Reuse one of the original cluster slots at each merge. This keeps the
    // matrix n × n rather than allocating a new row for every internal node.
    let mut matrix = vec![0.0; count * count];
    let mut index = 0;
    for row in 0..count {
        for column in row + 1..count {
            let distance = distances[index];
            matrix[row * count + column] = distance;
            matrix[column * count + row] = distance;
            index += 1;
        }
    }

    let mut active = vec![true; count];
    let mut sizes = vec![1_usize; count];
    let mut nodes = (0..count).map(Node::Leaf).collect::<Vec<_>>();
    let mut chain = Vec::with_capacity(count);
    let mut active_count = count;

    while active_count > 1 {
        if chain.is_empty() {
            chain.push(first_active(&active));
        }
        let current = *chain.last().expect("chain was initialized");
        let nearest = nearest_neighbour(current, &active, &matrix, count);

        if chain.len() >= 2 && nearest == chain[chain.len() - 2] {
            let one = chain.pop().expect("current chain entry");
            let two = chain.pop().expect("previous chain entry");
            merge_ward(
                one,
                two,
                &mut matrix,
                &mut active,
                &mut sizes,
                &mut nodes,
                count,
            );
            active_count -= 1;
        } else {
            chain.push(nearest);
        }
    }

    let root = first_active(&active);
    nodes.swap_remove(root)
}

fn first_active(active: &[bool]) -> usize {
    active
        .iter()
        .position(|is_active| *is_active)
        .expect("at least one cluster remains active")
}

fn nearest_neighbour(current: usize, active: &[bool], matrix: &[f64], count: usize) -> usize {
    active
        .iter()
        .enumerate()
        .filter(|(candidate, is_active)| *candidate != current && **is_active)
        .min_by(|(left, _), (right, _)| {
            matrix[current * count + *left]
                .total_cmp(&matrix[current * count + *right])
                .then_with(|| left.cmp(right))
        })
        .map(|(index, _)| index)
        .expect("an active cluster has a neighbour")
}

fn merge_ward(
    one: usize,
    two: usize,
    matrix: &mut [f64],
    active: &mut [bool],
    sizes: &mut [usize],
    nodes: &mut [Node],
    count: usize,
) {
    // Keeping the lower slot makes ties and resulting leaf order reproducible.
    let (keep, remove) = (one.min(two), one.max(two));
    let (one_size, two_size) = (sizes[keep] as f64, sizes[remove] as f64);
    let between = matrix[keep * count + remove];

    for other in 0..count {
        if other == keep || other == remove || !active[other] {
            continue;
        }
        let other_size = sizes[other] as f64;
        let keep_distance = matrix[keep * count + other];
        let remove_distance = matrix[remove * count + other];
        // Lance-Williams update for Ward's minimum-variance linkage.
        let squared = ((one_size + other_size) * keep_distance.powi(2)
            + (two_size + other_size) * remove_distance.powi(2)
            - other_size * between.powi(2))
            / (one_size + two_size + other_size);
        let distance = squared.max(0.0).sqrt();
        matrix[keep * count + other] = distance;
        matrix[other * count + keep] = distance;
    }

    let old_keep = std::mem::replace(&mut nodes[keep], Node::Leaf(keep));
    let old_remove = std::mem::replace(&mut nodes[remove], Node::Leaf(remove));
    nodes[keep] = Node::Merge(Box::new(old_keep), Box::new(old_remove));
    active[remove] = false;
    sizes[keep] += sizes[remove];
}

#[derive(Debug)]
enum Node {
    Leaf(usize),
    Merge(Box<Node>, Box<Node>),
}

#[derive(Clone, Copy)]
struct State {
    cost: f64,
    first: usize,
    last: usize,
    trace: usize,
}

#[derive(Clone, Copy)]
struct Join {
    left: (usize, bool),
    right: (usize, bool),
}

struct Candidate {
    cost: f64,
    first: usize,
    last: usize,
    join: Join,
}

fn bounded_leaf_order(node: &Node, distances: &[f64], count: usize, limit: usize) -> Vec<usize> {
    fn fill(
        node: &Node,
        distances: &[f64],
        count: usize,
        limit: usize,
        traces: &mut Vec<Join>,
    ) -> Vec<State> {
        let (left, right) = match node {
            Node::Leaf(i) => {
                return vec![State {
                    cost: 0.,
                    first: *i,
                    last: *i,
                    trace: *i,
                }];
            }
            Node::Merge(left_node, right_node) => (
                fill(left_node, distances, count, limit, traces),
                fill(right_node, distances, count, limit, traces),
            ),
        };
        let mut candidates = Vec::with_capacity(4 * left.len() * right.len());
        for left_order in left {
            for right_order in &right {
                for reverse_left in [false, true] {
                    for reverse_right in [false, true] {
                        let (left_start, left_end) = if reverse_left {
                            (left_order.last, left_order.first)
                        } else {
                            (left_order.first, left_order.last)
                        };
                        let (right_start, right_end) = if reverse_right {
                            (right_order.last, right_order.first)
                        } else {
                            (right_order.first, right_order.last)
                        };
                        let cost = left_order.cost
                            + right_order.cost
                            + condensed_distance(distances, count, left_end, right_start);
                        // A path and its reversal have equal cost. Store one orientation.
                        let (first, last, join) = if left_start < right_end {
                            (
                                left_start,
                                right_end,
                                Join {
                                    left: (left_order.trace, reverse_left),
                                    right: (right_order.trace, reverse_right),
                                },
                            )
                        } else {
                            (
                                right_end,
                                left_start,
                                Join {
                                    left: (right_order.trace, !reverse_right),
                                    right: (left_order.trace, !reverse_left),
                                },
                            )
                        };
                        candidates.push(Candidate {
                            cost,
                            first,
                            last,
                            join,
                        });
                    }
                }
            }
        }
        // Only the cheapest path for each endpoint pair can help its parent.
        candidates.sort_by(|left_candidate, right_candidate| {
            (left_candidate.first, left_candidate.last)
                .cmp(&(right_candidate.first, right_candidate.last))
                .then_with(|| left_candidate.cost.total_cmp(&right_candidate.cost))
        });
        candidates.dedup_by(|left_candidate, right_candidate| {
            left_candidate.first == right_candidate.first
                && left_candidate.last == right_candidate.last
        });
        candidates.sort_by(|left_candidate, right_candidate| {
            left_candidate
                .cost
                .total_cmp(&right_candidate.cost)
                .then_with(|| {
                    (left_candidate.first, left_candidate.last)
                        .cmp(&(right_candidate.first, right_candidate.last))
                })
        });
        candidates.truncate(limit);
        candidates
            .into_iter()
            .map(|candidate| {
                let trace = count + traces.len();
                traces.push(candidate.join);
                State {
                    cost: candidate.cost,
                    first: candidate.first,
                    last: candidate.last,
                    trace,
                }
            })
            .collect()
    }
    fn trace(index: usize, reverse: bool, count: usize, traces: &[Join], out: &mut Vec<usize>) {
        if index < count {
            out.push(index);
        } else {
            // Leaf indices occupy 0..count; later indices refer to stored joins.
            let Join {
                left: (left, reverse_left),
                right: (right, reverse_right),
            } = traces[index - count];
            if reverse {
                trace(right, !reverse_right, count, traces, out);
                trace(left, !reverse_left, count, traces, out);
            } else {
                trace(left, reverse_left, count, traces, out);
                trace(right, reverse_right, count, traces, out);
            }
        }
    }
    let mut traces = Vec::new();
    let states = fill(node, distances, count, limit, &mut traces);
    let mut result = Vec::with_capacity(count);
    trace(states[0].trace, false, count, &traces, &mut result);
    result
}

fn condensed_distance(distances: &[f64], count: usize, one: usize, two: usize) -> f64 {
    let (row, column) = (one.min(two), one.max(two));
    distances[row * (2 * count - row - 1) / 2 + column - row - 1]
}

#[cfg(test)]
mod tests {
    use super::ward_leaf_order;

    #[test]
    fn orders_two_close_pairs_next_to_each_other() {
        // (0, 1) and (2, 3) are the nearest pairs.
        let order = ward_leaf_order(&[1.0, 8.0, 8.0, 8.0, 8.0, 1.0], 4);
        assert_eq!(order.len(), 4);
        assert!(
            order
                .windows(2)
                .any(|pair| pair == [0, 1] || pair == [1, 0])
        );
        assert!(
            order
                .windows(2)
                .any(|pair| pair == [2, 3] || pair == [3, 2])
        );
    }

    #[test]
    fn orients_ward_branches_to_keep_linked_clusters_adjacent() {
        // Ward first forms (0, 1) and (2, 3). The old construction-order
        // traversal returned 0, 1, 2, 3, making the unlinked 1/2 pair
        // neighbours. Flipping the second branch preserves both groups and
        // makes its boundary the linked 1/3 pair instead.
        let order = ward_leaf_order(&[0.1, 1.0, 1.0, 1.0, 0.2, 0.1], 4);
        assert!(
            order
                .windows(2)
                .any(|pair| pair == [1, 3] || pair == [3, 1])
        );
        assert!(
            !order
                .windows(2)
                .any(|pair| pair == [1, 2] || pair == [2, 1])
        );
    }

    fn all_orders(node: &super::Node) -> Vec<Vec<usize>> {
        match node {
            super::Node::Leaf(i) => vec![vec![*i]],
            super::Node::Merge(left, right) => {
                let mut orders = Vec::new();
                for left in all_orders(left) {
                    for right in all_orders(right) {
                        let mut order = left.clone();
                        order.extend(&right);
                        orders.push(order.clone());
                        order.reverse();
                        orders.push(order);
                    }
                }
                orders
            }
        }
    }

    #[test]
    fn matches_exhaustive_search_when_unpruned_and_preserves_groups_when_pruned() {
        use super::Node::{Leaf, Merge};
        let pair = |a, b| Merge(Box::new(Leaf(a)), Box::new(Leaf(b)));
        let trees = [
            Merge(
                Box::new(Merge(Box::new(pair(0, 1)), Box::new(pair(2, 3)))),
                Box::new(Merge(Box::new(Leaf(4)), Box::new(pair(5, 6)))),
            ),
            (1..7).fold(Leaf(0), |tree, i| Merge(Box::new(tree), Box::new(Leaf(i)))),
        ];
        for tree in &trees {
            let possibilities = all_orders(tree);
            for seed in 0..12 {
                let distances = (0..21)
                    .map(|i| ((i * 17 + seed * 13 + i * seed) % 31) as f64 / 31.0)
                    .collect::<Vec<_>>();
                let cost = |order: &[usize]| {
                    order
                        .windows(2)
                        .map(|pair| super::condensed_distance(&distances, 7, pair[0], pair[1]))
                        .sum::<f64>()
                };
                let expected = possibilities
                    .iter()
                    .map(|order| cost(order))
                    .min_by(f64::total_cmp)
                    .unwrap();
                let actual = super::bounded_leaf_order(tree, &distances, 7, 32);
                assert!(possibilities.contains(&actual));
                for limit in [1, 4] {
                    let pruned = super::bounded_leaf_order(tree, &distances, 7, limit);
                    assert!(possibilities.contains(&pruned));
                    assert!(cost(&pruned) >= expected - 1e-12);
                }
                assert!((cost(&actual) - expected).abs() < 1e-12);
                assert_eq!(super::bounded_leaf_order(tree, &distances, 7, 32), actual);
            }
        }
        assert!(ward_leaf_order(&[], 0).is_empty());
        assert_eq!(ward_leaf_order(&[], 1), [0]);
    }
}
