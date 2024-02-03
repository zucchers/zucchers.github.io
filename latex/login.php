<?php
// Define your username and password
$username = "latex";
$password = "letmein";
if ($_POST['txtUsername'] != $username || $_POST['txtPassword'] != $password) {
?>

<title>Generatore di verifiche con soluzione</title>
<h1>Area protetta</h1>
Accesso consentito solo tramite credenziali<br>
<form name="form" method="post" action="<?php echo $_SERVER['PHP_SELF']; ?>">
    <p><label for="txtUsername">Nome utente:</label>
    <br /><input type="text" title="Enter your Username" name="txtUsername" /></p>
    <p><label for="txtpassword">Password:</label>
    <br /><input type="password" title="Enter your password" name="txtPassword" /></p>
    <p><input type="submit" name="Submit" value="Accedi" /></p>
</form>

<?php
}
else {
?>
<title>Generatore di verifiche</title>
<body>
<p>Gli esercizi ed i relativi codici sono scaricabili <a href="book.pdf#toolbar=0&amp;navpanes=0&amp;scrollbar=1&amp;page=1&amp;view=FitH" target="_blank">qui</a> (file pdf criptato e non stampabile, password necessaria per visualizzarlo)</p>

<form action="./processa.php" method="post">
<p><label for="txtTitolo">Titolo della verifica:</label><br />
<input type="text" title="Inserire il titolo della verifica" name="txtTitolo" /></p>
<p><label for="txtData">Data della verifica:</label><br />
<input type="text" title="Inserire la data della verifica" name="txtData" /></p>
<p><label for="txtClasse">Classe:</label><br />
<input type="text" title="Inserire la classe (es. 4F)" name="txtClasse" /></p>
<p><label for="txtLui">Lui:</label><br />
<input type="text" title="Inserire un nome maschile " name="txtLui" /></p>
<p><label for="txtLei">Lei:</label><br />
<input type="text" title="Inserire un nome femminile " name="txtLei" /></p>
<p><label for="txtLista">Lista dei codici degli esercizi (andare a caporiga):</label><br />
<TEXTAREA type="text" title="Inserire i codici andando a caporiga"  name="txtLista" rows="8" cols="23"></TEXTAREA><br>
<p></p>
<input type="submit" value="Genera verifica" />
</form>
</body>
</html>
<?php
}
?> 

