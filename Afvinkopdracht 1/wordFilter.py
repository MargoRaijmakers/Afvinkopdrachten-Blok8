# Margo Raijmakers
# 12-04-2021
# Afvinkopdracht 1: Basis Flask

from flask import Flask, render_template, request
import mysql.connector
import re

app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
def get_filter():
    """Deze functie zorgt ervoor dat er een zoekwoord opgegeven kan
    worden en dat de descriptions met het zoekwoord getoond worden.

    :return: het HTML bestand met de filter en data
    """
    if request.method == "POST":
        ftr = request.form.get("filter", "")
        conn = connector()
        ftr_data = use_filter(conn, ftr)
        return render_template("WebAppWordFilter.html", filter=ftr,
                               data=ftr_data)
    else:
        return render_template("WebAppWordFilter.html", filter="",
                               data="")


def connector():
    """Deze functie connect aan de database.

    :return conn: de connectie aan de database
    """
    conn = mysql.connector.connect(host="ensembldb.ensembl.org",
                                   user="anonymous",
                                   db="homo_sapiens_core_95_38")
    return conn


def use_filter(conn, ftr):
    """Deze functie filtert de data met het zoekwoord in de description.

    :param conn: de connectie aan de database
    :param ftr: de opgegeven filter
    :return ftr_data: de data met het zoekwoord in de description
    """
    cursor = conn.cursor()
    cursor.execute("select description from gene")
    rows = cursor.fetchall()
    ftr_data = []
    for row in rows:
        if re.search(ftr, str(row[0])):
            ftr_data.append(search_word_pos(str(row[0]), ftr))
    cursor.close()
    conn.close()
    return ftr_data


def search_word_pos(line, ftr):
    """Deze functie split de regel op het filterwoord.

    :param line: de regel
    :param ftr: de opgegeven filter
    :return: de gesplitte regel in een lijst
    """
    line_split = line.split(ftr, 1)
    line_split.insert(1, ftr)
    return line_split


if __name__ == "__main__":
    app.run()
