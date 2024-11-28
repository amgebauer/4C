#!/bin/bash
# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Exit the script at the first failure
set -e

find dependencies docker -not -wholename '*/trilinos_develop/*' -not -name 'README.md' -type f -exec sha1sum {} \; | sort | sha1sum | cut -c -8