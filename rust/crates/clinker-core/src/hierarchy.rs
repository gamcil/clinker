//! Agglomerative hierarchy construction for cluster ordering.
/// Return the leaves of a Ward-linkage dendrogram for a condensed distance
/// matrix. The nearest-neighbour-chain method uses O(n²) time and memory.
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
    let mut leaves = Vec::with_capacity(count);
    nodes[root].append_leaves(&mut leaves);
    leaves
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

impl Node {
    fn append_leaves(&self, leaves: &mut Vec<usize>) {
        match self {
            Self::Leaf(index) => leaves.push(*index),
            Self::Merge(left, right) => {
                left.append_leaves(leaves);
                right.append_leaves(leaves);
            }
        }
    }
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
}
