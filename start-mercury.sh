#!/bin/sh

echo "Add mercury web applications"
for filepath in webapps/*.ipynb; do
    echo "Add " $filepath
    poetry run mercury add $filepath
done

poetry run mercury run
