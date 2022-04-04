#!/bin/sh
ARCHIVE_DIR=$(head -n 1 ARCHIVE_DIR)
echo $ARCHIVE_DIR
cd $ARCHIVE_DIR
wget -r --no-clobber --no-parent --continue http://cdn.gea.esac.esa.int/Gaia/gedr3/gaia_source/
