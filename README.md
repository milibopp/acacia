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
- Arbitrary computational queries can be performed on the trees.
- Faster dimension-specific trees, e.g. an octree (not yet implemented)


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
