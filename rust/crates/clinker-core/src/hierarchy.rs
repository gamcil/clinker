//! Agglomerative hierarchy construction and leaf ordering for cluster ordering.

use std::collections::BTreeMap;

/// Return the leaves of a Ward-linkage dendrogram for a condensed distance
/// matrix.
///
/// The nearest-neighbour-chain construction uses O(n²) time and memory. Once
/// the tree is built, optimal leaf ordering chooses the orientation of each
/// branch that minimizes the sum of distances between displayed neighbours.
/// This retains Ward's groups while avoiding arbitrary branch boundaries that
/// can otherwise put two clusters with no links next to one another.
pub(crate) fn ward_leaf_order(distances: &[f64], count: usize) -> Vec<usize> {
    if count < 2 {
        return (0..count).collect();
    }
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
    optimal_leaf_order(&nodes[root], distances, count)
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

#[derive(Debug, Clone)]
struct OrderedLeaves {
    cost: f64,
    leaves: Vec<usize>,
}

/// Find the lowest-cost sequence for every possible pair of endpoints in a
/// subtree. The recurrence only joins a left and right child, so every result
/// remains a valid ordering of the Ward dendrogram.
fn optimal_leaf_order(node: &Node, distances: &[f64], count: usize) -> Vec<usize> {
    let states = optimal_leaf_order_states(node, distances, count);
    states
        .into_values()
        .min_by(|left, right| compare_orders(left, right))
        .expect("a dendrogram contains a leaf")
        .leaves
}

fn optimal_leaf_order_states(
    node: &Node,
    distances: &[f64],
    count: usize,
) -> BTreeMap<(usize, usize), OrderedLeaves> {
    match node {
        Node::Leaf(index) => BTreeMap::from([(
            (*index, *index),
            OrderedLeaves {
                cost: 0.0,
                leaves: vec![*index],
            },
        )]),
        Node::Merge(left, right) => {
            let left_states = optimal_leaf_order_states(left, distances, count);
            let right_states = optimal_leaf_order_states(right, distances, count);
            let mut states = BTreeMap::new();

            for (&(left_start, left_end), left_order) in &left_states {
                for (&(right_start, right_end), right_order) in &right_states {
                    let cost = left_order.cost
                        + condensed_distance(distances, count, left_end, right_start)
                        + right_order.cost;
                    let mut leaves = left_order.leaves.clone();
                    leaves.extend_from_slice(&right_order.leaves);
                    keep_best(
                        &mut states,
                        (left_start, right_end),
                        OrderedLeaves { cost, leaves },
                    );

                    let cost = right_order.cost
                        + condensed_distance(distances, count, right_end, left_start)
                        + left_order.cost;
                    let mut leaves = right_order.leaves.clone();
                    leaves.extend_from_slice(&left_order.leaves);
                    keep_best(
                        &mut states,
                        (right_start, left_end),
                        OrderedLeaves { cost, leaves },
                    );
                }
            }
            states
        }
    }
}

fn condensed_distance(distances: &[f64], count: usize, one: usize, two: usize) -> f64 {
    if one == two {
        return 0.0;
    }
    let (row, column) = (one.min(two), one.max(two));
    let index = row * (2 * count - row - 1) / 2 + (column - row - 1);
    distances[index]
}

fn keep_best(
    states: &mut BTreeMap<(usize, usize), OrderedLeaves>,
    endpoints: (usize, usize),
    candidate: OrderedLeaves,
) {
    match states.get(&endpoints) {
        Some(current) if compare_orders(current, &candidate).is_le() => {}
        _ => {
            states.insert(endpoints, candidate);
        }
    }
}

fn compare_orders(left: &OrderedLeaves, right: &OrderedLeaves) -> std::cmp::Ordering {
    left.cost
        .total_cmp(&right.cost)
        .then_with(|| left.leaves.cmp(&right.leaves))
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
}
