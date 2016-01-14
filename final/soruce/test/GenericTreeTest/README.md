## README.md ##

1. `make`
2. `make testdata`
3. `make test`
4. `make clean`

### Result ###

| sizeof(LinearTreeNode) bytes \ Traversal | Recursive | Loop |
|--------|--------|--------|
|  32    |  6.049s|  5.628s|
|  44    |  6.651s|  6.817s|
|  60    |  7.460s|  6.888s|
|  92    |  9.361s|  9.271s|
| 156    | 16.844s| 16.694s|
| 220    | 25.294s| 27.031s|
| 284    | 28.181s| 30.900s|
| 540    | 28.560s| 33.707s|