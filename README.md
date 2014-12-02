`acacia` is my attempt at writing a generic spatial tree library including
binary trees, quadtrees, octrees and so on. The idea is to associate a generic
data field to each node to be computed using a suitable closure. One can then
pass arbitrary computational queries to the tree.

This also allows reusing the same tree structure for different queries without
mixing their implementation within the same type.

The current state of the project is very experimental. It works, as the tests
indicate, but I intend to change both the API and the internals significantly
to improve on design and performance. There is also no documentation at the
moment.


# License

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
