#!/bin/bash
echo "Avvio il sito con il path corretto"
export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python3.10/site-packages/
echo "Sito avviato"
python3.10 app.py