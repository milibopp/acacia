[![Build Status](https://img.shields.io/travis/aepsil0n/acacia.svg)](https://travis-ci.org/aepsil0n/acacia)
[![](https://img.shields.io/crates/v/acacia.svg)](https://crates.io/crates/acacia)

*acacia* is a spatial tree library written in Rust. It is generic over the
dimension of the partitioned space and thus supports binary trees, quadtrees,
octrees, etc. The intended goal is to implement these features as fast and
covering as many use cases as possible without sacrificing abstraction.

The current state of the project is very experimental. It works and has ample
test coverage, but both the API and the internals will likely change in the
future to improve interface and performance.


## Features

- Tree construction from a simple iterator.
- Associate data to a tree during construction using closures.
- Perform arbitrary computational queries on the trees.
- Faster dimension-specific trees, e.g. an octree (not yet implemented)


## Documentation

The documentation is hosted on [Rust CI][rustci_acacia] and should be up-to-date
with the master branch.


## Example: N-body gravity calculation

Gravity calculations are a fairly common example for speeding up calculations
with spatial trees. This is a simple example to calculate the gravitational
acceleration at a given point between a set of gravitating particles. The code
presented here is an excerpt from a complete example you can find in the
directory `example/gravity`.

The tree can be constructed from an iterator over particles and some data about
its geometry.

```rust
let tree = Tree::new(
    particles.iter(),
    Ncube::new(origin, 11.0),
```

Note, that the particles implement the `Position` trait used to define the
notion that a type has a position.

Next we need a couple of values and closures to associate data with each node in
the tree: a value for empty nodes, a closure that assigns a value to a leaf node
given the object stored in it, and a closure that combines two pieces of
associated data to compute values for branch nodes.

```rust
    (origin, 0.0),
    &|obj| (obj.position, obj.mass),
    &|&(com1, m1), &(com2, m2)|
        if m1 + m2 > 0.0 {(
            origin + (com1.to_vec() * m1 + com2.to_vec() * m2) / (m1 + m2),
            m1 + m2,
        )}
        else {
            (origin, 0.0)
        }
);
```

The associated data in this example is a tuple made from the center of mass and
the total mass of a node.

Now we can issue a computational query to the tree by passing in two more
closures to its `query_data` method: the first one serves as a criterion for
recursion. If a branch node passes this, the query continues on its children.
The second one collects force terms from each piece of associated data the tree
encounters during this recursion.

```rust
let mut tree_gravity: Vec3<f64> = zero();
tree.query_data(
    &|node| {
        let &(ref center_of_mass, _) = node.data();
        let d = test_point.dist(center_of_mass);
        let delta: f64 = node.partition().center().dist(center_of_mass);
        d < 2.0 * node.partition().width() + delta
    },

    &mut |&(center_of_mass, mass)| {
        tree_gravity = tree_gravity + newton(mass, center_of_mass, test_point);
    },
);
```


## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/.


[rustci_acacia]: http://www.rust-ci.org/aepsil0n/acacia/doc/acacia/
